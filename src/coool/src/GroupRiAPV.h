#ifndef __GroupRiAPV__
#define __GroupRiAPV__

#include "Group.h"
#include <map>


class GroupRiAPV : public Group {

 private:

  typedef struct { double sum, sum2; int nb; } stat_s ;
  class PdRefType
                   { public: int xpos, ypos, plindex;
                     PdRefType() {}
                     PdRefType(int x, int y, int pl): xpos(x), ypos(y), plindex(pl) {}
		     PdRefType(const PdRefType& rf): xpos(rf.xpos), ypos(rf.ypos), plindex(rf.plindex) {}
                   };

  /// number of planes in the group
  int fNbPlanes;

  /// number of planes used by the group
  int fNbPlanesUsed;

  TH2F *fHch2Dall, *fHchamp2Dall;
  TH2D *fHavgamp2Dall, *fHsigamp2Dall;
  TH2D *fHavgampCM2Dall, *fHsigampCM2Dall;
  TH2D *fHevtamp2Dall, *fHevtampcut2Dall;
  TH2F *fHapvvsa2meanall;

  /// map to store characteristics of each PD
  std::map<unsigned int, PdRefType> pdRef;

  /// store a new PD in pdRef
  void addPdRef(unsigned int i, const PdRefType& rf);


 public:
  GroupRiAPV(const char* name);

  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif


  ClassDef(GroupRiAPV,0)
};


#endif
