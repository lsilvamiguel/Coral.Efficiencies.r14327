#include "CsTimingInfo.h"

#include <stdio.h>
#include <cassert>
#include <cstdio>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "CsTimingInfoFormat.h"
#include "CsEvent.h"
#include "CsOpt.h"
#include "CsErrLog.h"


// initialization of singleton pointer
CsTimingInfo* CsTimingInfo::instance = NULL;


CsTimingInfo::CsTimingInfo() : fMode(1),fValidFrom(-999),fValidTo(-999),fRun(0) {
  // This chould never be called twice
  assert( !instance );
  for( unsigned int spill=0; spill<300; ++spill) {
    unsigned int checksum(0);
    for(int srcid=0; srcid<6; ++srcid) {
      for(int port=0; port<8; ++port) {
        for(int ch_grp=0; ch_grp<8; ++ch_grp) {
          fArray[spill][srcid][port][ch_grp]=0;     // no time offset
          fEvnr[spill][srcid][port][ch_grp]=9e6;     //number higher than maximum number of events per spill,
          //in spill time corrections will be skipped by that
        }
      }
    }
  }
  if(CsEvent::Instance()->getRunNumber()<=82363 && CsEvent::Instance()->getRunNumber()>=72003)
    fType=0;
  else
    fType=1;
};

CsTimingInfo::~CsTimingInfo() {};

CsTimingInfo* CsTimingInfo::Instance() {
  if( !instance ) instance = new CsTimingInfo();
  return instance;
};

float CsTimingInfo::GetEC02Shift(const unsigned int srcid, const unsigned int port, const unsigned int chip, const unsigned int channel) {
  const static unsigned int min_srcid(616);

  // Check that parameters are in range
  if( srcid < min_srcid ||
      srcid > min_srcid+5 ||
      port > 7 ||
      chip > 3 ||
      channel > 16){
    return 0.;
  }
  unsigned int chgrp = (chip<<1) + (channel>>3);

  //Check if data is valid
  if(CsEvent::Instance()->getBurstNumber()-1<300)
    if(CsEvent::Instance()->getEventNumberInBurst()>=fEvnr[CsEvent::Instance()->getBurstNumber()-1][srcid-min_srcid][port][chgrp])
      return fArray[CsEvent::Instance()->getBurstNumber()-1][srcid-min_srcid][port][chgrp]*12.86;       //shift time by N clock cycles of 80 MHz

  return 0;

};

int CsTimingInfo::ReadEC02ShiftFile(std::string s) {

  int runNmb=CsEvent::Instance()->getRunNumber();

  std::istringstream is(s);
  std::string str;
  while( getline(is,str) ) {
    if(str=="" || str.c_str()[0]=='#') continue;
    int run, spill, location, evnr;
    int srcid,port,ch_grp;
    int dir;
    int ret = sscanf(str.c_str(),"%d %d %d %d %d",&run,&spill,&location,&evnr,&dir);
    assert( ret == 5 );
    if(run == runNmb || run==0) {
      location--;
      srcid=location>>6;
      port=((location&0x3f)>>3);
      ch_grp=location&7;
      if(run!=0) {
        fArray[spill][srcid][port][ch_grp]=-dir;
        fEvnr[spill][srcid][port][ch_grp]=evnr;
      }
      else
        for(int h=0;h<300;h++) {
          fArray[h][srcid][port][ch_grp]=-dir;
          fEvnr[h][srcid][port][ch_grp]=evnr;
        }
      jump_evnr[spill].push_back(evnr);
    };


  }

  return 0;
};

unsigned int CsTimingInfo::GetDistToJump() {
  //returns the distance to the closest in spill time jump
  //it gives the possibility to cut on this value
  int dist=(1<<25);
  int spill=CsEvent::Instance()->getBurstNumber()-1;
  if(spill<300)
    for(unsigned int i=0;i<jump_evnr[spill].size();i++) {
      dist = std::min(dist,std::abs((int)CsEvent::Instance()->getEventNumberInBurst()-jump_evnr[spill][i]));
    }
  return dist;
};

unsigned int CsTimingInfo::GetTriggerAffected() {
  //only for 2009 Primakoff production
  static TriggerAddr* TriggerMap;
  static bool spill_checked[200]={false};
  static bool first(true);
  static unsigned int affectedMask[200]={0};
  if(first) {
    first=false;
    TriggerMap =new TriggerAddr[24];
    TriggerMap[0].srcid=616;  TriggerMap[0].port=4; TriggerMap[0].chip=2;  TriggerMap[0].channel=0;
    TriggerMap[1].srcid=616;  TriggerMap[1].port=4; TriggerMap[1].chip=3;  TriggerMap[1].channel=1;
    TriggerMap[2].srcid=616;  TriggerMap[2].port=5; TriggerMap[2].chip=0;  TriggerMap[2].channel=0;
    TriggerMap[3].srcid=616;  TriggerMap[3].port=5; TriggerMap[3].chip=1;  TriggerMap[3].channel=1;
    TriggerMap[4].srcid=616;  TriggerMap[4].port=5; TriggerMap[4].chip=2;  TriggerMap[4].channel=0;
    TriggerMap[5].srcid=616;  TriggerMap[5].port=5; TriggerMap[5].chip=3;  TriggerMap[5].channel=1;
    TriggerMap[6].srcid=616;  TriggerMap[6].port=6; TriggerMap[6].chip=0;  TriggerMap[6].channel=0;
    TriggerMap[7].srcid=616;  TriggerMap[7].port=6; TriggerMap[7].chip=1;  TriggerMap[7].channel=1;
    TriggerMap[8].srcid=616;  TriggerMap[8].port=6; TriggerMap[8].chip=2;  TriggerMap[8].channel=0;
    TriggerMap[9].srcid=616;  TriggerMap[9].port=6; TriggerMap[9].chip=3;  TriggerMap[9].channel=1;
    TriggerMap[10].srcid=618; TriggerMap[10].port=1;  TriggerMap[10].chip=0; TriggerMap[10].channel=0;
    TriggerMap[11].srcid=618; TriggerMap[11].port=1;  TriggerMap[11].chip=1; TriggerMap[11].channel=1;
    TriggerMap[12].srcid=618; TriggerMap[12].port=1;  TriggerMap[12].chip=2; TriggerMap[12].channel=0;
    TriggerMap[13].srcid=618; TriggerMap[13].port=1;  TriggerMap[13].chip=3; TriggerMap[13].channel=1;
    TriggerMap[14].srcid=618; TriggerMap[14].port=0;  TriggerMap[14].chip=0; TriggerMap[14].channel=0;
    TriggerMap[15].srcid=618; TriggerMap[15].port=0;  TriggerMap[15].chip=1; TriggerMap[15].channel=1;
    TriggerMap[16].srcid=618; TriggerMap[16].port=0;  TriggerMap[16].chip=2; TriggerMap[16].channel=0;
    TriggerMap[17].srcid=618; TriggerMap[17].port=0;  TriggerMap[17].chip=3; TriggerMap[17].channel=1;
    TriggerMap[18].srcid=621; TriggerMap[18].port=6;  TriggerMap[18].chip=0; TriggerMap[18].channel=0;
    TriggerMap[19].srcid=621; TriggerMap[19].port=6;  TriggerMap[19].chip=1; TriggerMap[19].channel=1;
    TriggerMap[20].srcid=621; TriggerMap[20].port=6;  TriggerMap[20].chip=2; TriggerMap[20].channel=0;
    TriggerMap[21].srcid=621; TriggerMap[21].port=6;  TriggerMap[21].chip=3; TriggerMap[21].channel=1;
    TriggerMap[22].srcid=621; TriggerMap[22].port=5;  TriggerMap[22].chip=0; TriggerMap[22].channel=0;
    TriggerMap[23].srcid=621; TriggerMap[23].port=5;  TriggerMap[23].chip=1; TriggerMap[23].channel=1;
  }
  unsigned int spill=CsEvent::Instance()->getBurstNumber()-1;
  if(spill<200) {
    if(!spill_checked[spill]) {
      spill_checked[spill]=true;
      affectedMask[spill]=0;
      for(unsigned int i=0;i<24;i++) {
        float shift=GetEC02Shift(TriggerMap[i].srcid, TriggerMap[i].port, TriggerMap[i].chip, TriggerMap[i].channel<<3);
        if(shift!=0) {
          affectedMask[spill]|=(1<<i);
        }
      }
    }
    return affectedMask[spill];
  }
  else
    return 0;
};

  unsigned int CsTimingInfo::GetTimingInfo() {
    if(fType==0)
      return GetTriggerAffected();
    if(fType==1)
      return GetDistToJump();
    else
      return (1<<25);     //could be masked out for Primakoff data,
    //for in spill timing corrections it's representing infinity distance to the next jump
  }
