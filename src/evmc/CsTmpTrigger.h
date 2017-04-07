#ifndef CsTmpTrigger_h
#define CsTmpTrigger_h

#include <list>
#include "CsGeant3.h"
#include "CsMCParticle.h"
#include "CsMCTrack.h"

typedef struct { 
  float qq2;
  float y;
  float phi;
  int pp;
  int l;
  int p;
  int unp;    
} TrigVar;

class CsTmpTrigger{

 public:
 
//methods
  bool CheckTmpTrigger();  //!< checks trigger  
  bool ReadTriggerMap();   //!< read kinematical trigger map

 private:
  std::list<TrigVar> tryg;    //!< maplist
   int conf;             //!< configuration=1 or 2001  
   unsigned int noTrig;  //!< number of not triggered events

};

#endif //CsTmpTrigger_h
