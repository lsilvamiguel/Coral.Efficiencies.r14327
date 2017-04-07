//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Sat Sep 15 00:59:48 2001 by ROOT version3.00/06)
//   from TTree T/COMPASS monitoring
//   found on file: ../root/mm.root
//////////////////////////////////////////////////////////


#ifndef MMan_h
#define MMan_h

#include <string>
#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMinuit.h"

#include "TrackMaker.h"
#include "Tracker.h"
#include "Doublet.h"
#include "Station.h"

#define MAX_MULT 1024*8


class MMan {

  private :
    bool fHitsInTree;
    int fMult;
    bool fAllocated;
    float fChi2;
    bool fDoTree;
    std::string fDetdat;
    TMinuit *fMinuit;
    int fNpar;
    std::vector<bool> fEventMask;
    int fNeventOpt;

    std::vector<TH1*> fHresSave;
    std::vector<TH1*> fHproSave;

  public :
    
    TrackMaker *fTracker;

   TFile          *fFile;
   TTree         *fChain;   //pointer to the analyzed TTree or TChain
   TTree          *fOutTree;

   Int_t           fCurrent; //current Tree number in a TChain
//Declaration of leaves types
   int           runno;
   int           spill;
   int           evtinrun;
   int           evtinspl;
   int           timesec;
   int           timemic;
   int             TT_hits;
   float           TT_ch[MAX_MULT];
   float           TT_t[MAX_MULT];



   // MICROMEGAS

   int             MM03Y1___hhits;
   float           MM03Y1___hch[MAX_MULT];
   float           MM03Y1___ht[MAX_MULT];
   float           MM03Y1___htot[MAX_MULT];
   float           MM03Y1___hs[MAX_MULT];
   float           MM03Y1___hres[MAX_MULT];

   int             MM03X1___hhits;
   float           MM03X1___hch[MAX_MULT];
   float           MM03X1___ht[MAX_MULT];
   float           MM03X1___htot[MAX_MULT];
   float           MM03X1___hs[MAX_MULT];
   float           MM03X1___hres[MAX_MULT];

   int             MM03U1___hhits;
   float           MM03U1___hch[MAX_MULT];
   float           MM03U1___ht[MAX_MULT];
   float           MM03U1___htot[MAX_MULT];
   float           MM03U1___hs[MAX_MULT];
   float           MM03U1___hres[MAX_MULT];

   int             MM03V1___hhits;
   float           MM03V1___hch[MAX_MULT];
   float           MM03V1___ht[MAX_MULT];
   float           MM03V1___htot[MAX_MULT];
   float           MM03V1___hs[MAX_MULT];
   float           MM03V1___hres[MAX_MULT];

   int             MM02Y1___hhits;
   float           MM02Y1___hch[MAX_MULT];
   float           MM02Y1___ht[MAX_MULT];
   float           MM02Y1___htot[MAX_MULT];
   float           MM02Y1___hs[MAX_MULT];
   float           MM02Y1___hres[MAX_MULT];

   int             MM02X1___hhits;
   float           MM02X1___hch[MAX_MULT];
   float           MM02X1___ht[MAX_MULT];
   float           MM02X1___htot[MAX_MULT];
   float           MM02X1___hs[MAX_MULT];
   float           MM02X1___hres[MAX_MULT];

   int             MM02U1___hhits;
   float           MM02U1___hch[MAX_MULT];
   float           MM02U1___ht[MAX_MULT];
   float           MM02U1___htot[MAX_MULT];
   float           MM02U1___hs[MAX_MULT];
   float           MM02U1___hres[MAX_MULT];

   int             MM02V1___hhits;
   float           MM02V1___hch[MAX_MULT];
   float           MM02V1___ht[MAX_MULT];
   float           MM02V1___htot[MAX_MULT];
   float           MM02V1___hs[MAX_MULT];
   float           MM02V1___hres[MAX_MULT];

   int             MM01Y1___hhits;
   float           MM01Y1___hch[MAX_MULT];
   float           MM01Y1___ht[MAX_MULT];
   float           MM01Y1___htot[MAX_MULT];
   float           MM01Y1___hs[MAX_MULT];
   float           MM01Y1___hres[MAX_MULT];

   int             MM01X1___hhits;
   float           MM01X1___hch[MAX_MULT];
   float           MM01X1___ht[MAX_MULT];
   float           MM01X1___htot[MAX_MULT];
   float           MM01X1___hs[MAX_MULT];
   float           MM01X1___hres[MAX_MULT];

   int             MM01U1___hhits;
   float           MM01U1___hch[MAX_MULT];
   float           MM01U1___ht[MAX_MULT];
   float           MM01U1___htot[MAX_MULT];
   float           MM01U1___hs[MAX_MULT];
   float           MM01U1___hres[MAX_MULT];

   int             MM01V1___hhits;
   float           MM01V1___hch[MAX_MULT];
   float           MM01V1___ht[MAX_MULT];
   float           MM01V1___htot[MAX_MULT];
   float           MM01V1___hs[MAX_MULT];
   float           MM01V1___hres[MAX_MULT];

//List of branches
   TBranch        *b_runno;
   TBranch        *b_spill;
   TBranch        *b_evtinrun;
   TBranch        *b_evtinspl;
   TBranch        *b_timesec;
   TBranch        *b_timemic;
   TBranch        *b_TT_hits;
   TBranch        *b_TT_ch;
   TBranch        *b_TT_t;

   TBranch        *b_MM03Y1___hhits;
   TBranch        *b_MM03Y1___hch;
   TBranch        *b_MM03Y1___ht;
   TBranch        *b_MM03Y1___htot;
   TBranch        *b_MM03Y1___hs;
   TBranch        *b_MM03Y1___hres;

   TBranch        *b_MM03X1___hhits;
   TBranch        *b_MM03X1___hch;
   TBranch        *b_MM03X1___ht;
   TBranch        *b_MM03X1___htot;
   TBranch        *b_MM03X1___hs;
   TBranch        *b_MM03X1___hres;

   TBranch        *b_MM03U1___hhits;
   TBranch        *b_MM03U1___hch;
   TBranch        *b_MM03U1___ht;
   TBranch        *b_MM03U1___htot;
   TBranch        *b_MM03U1___hs;
   TBranch        *b_MM03U1___hres;

   TBranch        *b_MM03V1___hhits;
   TBranch        *b_MM03V1___hch;
   TBranch        *b_MM03V1___ht;
   TBranch        *b_MM03V1___htot;
   TBranch        *b_MM03V1___hs;
   TBranch        *b_MM03V1___hres;

   TBranch        *b_MM02Y1___hhits;
   TBranch        *b_MM02Y1___hch;
   TBranch        *b_MM02Y1___ht;
   TBranch        *b_MM02Y1___htot;
   TBranch        *b_MM02Y1___hs;
   TBranch        *b_MM02Y1___hres;

   TBranch        *b_MM02X1___hhits;
   TBranch        *b_MM02X1___hch;
   TBranch        *b_MM02X1___ht;
   TBranch        *b_MM02X1___htot;
   TBranch        *b_MM02X1___hs;
   TBranch        *b_MM02X1___hres;

   TBranch        *b_MM02U1___hhits;
   TBranch        *b_MM02U1___hch;
   TBranch        *b_MM02U1___ht;
   TBranch        *b_MM02U1___htot;
   TBranch        *b_MM02U1___hs;
   TBranch        *b_MM02U1___hres;

   TBranch        *b_MM02V1___hhits;
   TBranch        *b_MM02V1___hch;
   TBranch        *b_MM02V1___ht;
   TBranch        *b_MM02V1___htot;
   TBranch        *b_MM02V1___hs;
   TBranch        *b_MM02V1___hres;

   TBranch        *b_MM01Y1___hhits;
   TBranch        *b_MM01Y1___hch;
   TBranch        *b_MM01Y1___ht;
   TBranch        *b_MM01Y1___htot;
   TBranch        *b_MM01Y1___hs;
   TBranch        *b_MM01Y1___hres;

   TBranch        *b_MM01X1___hhits;
   TBranch        *b_MM01X1___hch;
   TBranch        *b_MM01X1___ht;
   TBranch        *b_MM01X1___htot;
   TBranch        *b_MM01X1___hs;
   TBranch        *b_MM01X1___hres;

   TBranch        *b_MM01U1___hhits;
   TBranch        *b_MM01U1___hch;
   TBranch        *b_MM01U1___ht;
   TBranch        *b_MM01U1___htot;
   TBranch        *b_MM01U1___hs;
   TBranch        *b_MM01U1___hres;

   TBranch        *b_MM01V1___hhits;
   TBranch        *b_MM01V1___hch;
   TBranch        *b_MM01V1___ht;
   TBranch        *b_MM01V1___htot;
   TBranch        *b_MM01V1___hs;
   TBranch        *b_MM01V1___hres;



   MMan(TFile *file, TTree *tree, bool fHitsInTree=true, int mult=0, float chi2=0);

   // will chain all chunks for this run 
   MMan(const char* dir,int run, const char* detdat, bool fHitsInTree=true, int mult=0, float chi2=0);

   ~MMan();

   void Reset();

   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop(int nentries, bool select);
   Bool_t Notify();
   void   Show(Int_t entry = -1);

   void Allocate(int run);
   bool TreatEvent();
   void Entry(Int_t i);
   void Test();

   void Efficiency(const char* sel,const char* outfile,int nevent=0,const char* opt=""); 
   void PlotTracking(int n);
   TH1F* PlotResiduals(int id,int n,const char* opt="");
   void Center(int id,int n);

   void EffVsChi2(const char* outfile, float chi2min, float chi2max, float ste);

   void TreeOn(bool dotree) {fDoTree = dotree;}
   
   // find time cuts for all planes (n sigmas)
   void TimeCuts(int n);

   // geometry optimization
   void InitMinuit(int nevent, float dpos); 

   void Fcn(int &npar, double *gin, double &f, double *par, int iflag);
   
   ClassDef(MMan,0)
};



#endif









