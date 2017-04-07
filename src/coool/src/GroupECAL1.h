#ifndef __GroupECAL1__
#define __GroupECAL1__

#include "Group.h"
#include "PlaneFEM.h"
#include "PlaneECAL1.h"
#include <vector>

using namespace std;

class GroupECAL1 : public Group {

 private:

  typedef struct { double sum, sum2; int nb; } stat_s ;
  class PdRefType
                   { public: int xpos, ypos, plindex;
                     PdRefType() {}
                     PdRefType(int x, int y, int pl): xpos(x), ypos(y), plindex(pl) {}
		     PdRefType(const PdRefType& rf): xpos(rf.xpos), ypos(rf.ypos), plindex(rf.plindex) {}
                   };

  /// FEM and ECAL planes present
  
  bool fPlanesOK;
  double RefLedMin,RefLedMax;

  const PlaneECAL1* fPlaneECAL1;
  const PlaneFEM* fPlaneFEM;

  TH2F *fHch2Dall, *fHchamp2Dall;
  TH2D *fHavgamp2Dall, *fHsigamp2Dall;
  TH2D *fHavgampCM2Dall, *fHsigampCM2Dall;
  TH2D *fHevtamp2Dall, *fHevtampcut2Dall;
  TH2F *fHapvvsa2meanall;
  TH1F *hmnled;
  TProfile *hLedProf, *href;
  vector<TH1F*> hLasAmplrenF;

 public:
  GroupECAL1(const char* name);

  void Init();

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif


  ClassDef(GroupECAL1,0)
};


#endif
