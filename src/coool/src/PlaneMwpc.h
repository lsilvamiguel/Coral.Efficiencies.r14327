#ifndef __PlaneMwpc__
#define __PlaneMwpc__

#include "Plane1V.h"

#include "TTree.h"
#include "TFile.h"

class Monitor;

/// Plane for MWPC's

class PlaneMwpc : public Plane1V {

 int fNstation;
 int fNplane;
 int fLowTime;
 int fHiTime;
 int fNwires;

  #ifndef __CINT__
  void StoreDigit(CS::Chip::Digit* digit);
  #endif

 private:
  //int nevent;
  //std::vector<cluster> clsv;

  /// Variables
  //Variable *fVcls;

  /// map of digits. indexed by channel number (ie ordered !). Replace with a set ??
  std::map<int,CDigit1> fDigits;

  /// variables 
  Variable *fVcch, *fVct, *fVcs;

  /// Histograms
  TH1F *fHchits,*fHcch,*fHct,*fHcs, *fHns;

  TProfile* fPeff;

  /// number of clusters passing the cuts
  int fNclustKept;

 public:
  
  PlaneMwpc(const char *detname, int nchan, int center, int width);
  
  ~PlaneMwpc() {}

//   void StoreDigit(int channel, int time);

#ifndef __CINT__
  void EndEvent(const CS::DaqEvent &event);
#endif
  /// book histograms and branch the tree
  void Init(TTree* tree =0);

  /// reset before each event
  void Reset();

  /// resets qll histograms
  void ResetHistograms();

  /// do the clustering 
  void Clusterize();
  
  ClassDef(PlaneMwpc,0)
};

#endif




