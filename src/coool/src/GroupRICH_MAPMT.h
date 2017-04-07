#ifndef __GroupRICH_MAPMT__
#define __GroupRICH_MAPMT__
#include "Group.h"
#include "TStyle.h"

/* ----------------------------------------------------------------------
   Maintainer & author: Frank Nerling 
   New Class for monitoring  RICH_MAPMT: (x, y, time), created April 2006
   Here: Group of all 4 MAPMT Quadrants  u
   Updated Version: 30.08.2006
   ---------------------------------------------------------------------- */

class GroupRICH_MAPMT : public Group {

 private:

  int fNbPlanes;

  TH2F *fHch2D, *fHch2D_cc, *fHch2D_cc_b,*fHchamp2D;
  TH2D *fHrc1, *fHsigamp2D, *fHavgampCM2D, *fHsigampCM2D;
  TH2D *fHevtampCM2D, *fHevtampCMcut2D;

 public:
  GroupRICH_MAPMT(const char* name): Group(name) {}

  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif


  ClassDef(GroupRICH_MAPMT,0)
};


#endif
