//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Sat Sep 15 00:59:48 2001 by ROOT version3.00/06)
//   from TTree T/COMPASS monitoring
//   found on file: ../root/mm.root
//////////////////////////////////////////////////////////


#ifndef MMan_h
#define MMan_h

#include <string>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

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
    string fDetdat;

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

   // MWPCS 
   int             PB01U1___chits;
   float           PB01U1___cch[MAX_MULT];
   float           PB01U1___ct[MAX_MULT];

   int             PB01X1___chits;
   float           PB01X1___cch[MAX_MULT];
   float           PB01X1___ct[MAX_MULT];

   int             PB02V1___chits;
   float           PB02V1___cch[MAX_MULT];
   float           PB02V1___ct[MAX_MULT];

   int             PB03U1___chits;
   float           PB03U1___cch[MAX_MULT];
   float           PB03U1___ct[MAX_MULT];

   int             PB03X1___chits;
   float           PB03X1___cch[MAX_MULT];
   float           PB03X1___ct[MAX_MULT];

   int             PB04V1___chits;
   float           PB04V1___cch[MAX_MULT];
   float           PB04V1___ct[MAX_MULT];

   int             PB05U1___chits;
   float           PB05U1___cch[MAX_MULT];
   float           PB05U1___ct[MAX_MULT];

   int             PB05X1___chits;
   float           PB05X1___cch[MAX_MULT];
   float           PB05X1___ct[MAX_MULT];

   int             PB06V1___chits;
   float           PB06V1___cch[MAX_MULT];
   float           PB06V1___ct[MAX_MULT];


   // MICROMEGAS

/*     int             MM03Y1___hhits; */
/*     float           MM03Y1___hch[MAX_MULT]; */
/*     float           MM03Y1___ht[MAX_MULT]; */
/*     float           MM03Y1___htot[MAX_MULT]; */
/*     float           MM03Y1___hs[MAX_MULT]; */

/*     int             MM03X1___hhits; */
/*     float           MM03X1___hch[MAX_MULT]; */
/*     float           MM03X1___ht[MAX_MULT]; */
/*     float           MM03X1___htot[MAX_MULT]; */
/*     float           MM03X1___hs[MAX_MULT]; */

/*     int             MM03U1___hhits; */
/*     float           MM03U1___hch[MAX_MULT]; */
/*     float           MM03U1___ht[MAX_MULT]; */
/*     float           MM03U1___htot[MAX_MULT]; */
/*     float           MM03U1___hs[MAX_MULT]; */

/*     int             MM03V1___hhits; */
/*     float           MM03V1___hch[MAX_MULT]; */
/*     float           MM03V1___ht[MAX_MULT]; */
/*     float           MM03V1___htot[MAX_MULT]; */
/*     float           MM03V1___hs[MAX_MULT]; */

/*     int             MM01Y1___hhits; */
/*     float           MM01Y1___hch[MAX_MULT]; */
/*     float           MM01Y1___ht[MAX_MULT]; */
/*     float           MM01Y1___htot[MAX_MULT]; */
/*     float           MM01Y1___hs[MAX_MULT]; */

/*     int             MM01X1___hhits; */
/*     float           MM01X1___hch[MAX_MULT]; */
/*     float           MM01X1___ht[MAX_MULT]; */
/*     float           MM01X1___htot[MAX_MULT]; */
/*     float           MM01X1___hs[MAX_MULT]; */

/*     int             MM01U1___hhits; */
/*     float           MM01U1___hch[MAX_MULT]; */
/*     float           MM01U1___ht[MAX_MULT]; */
/*     float           MM01U1___htot[MAX_MULT]; */
/*     float           MM01U1___hs[MAX_MULT]; */

/*     int             MM01V1___hhits; */
/*     float           MM01V1___hch[MAX_MULT]; */
/*     float           MM01V1___ht[MAX_MULT]; */
/*     float           MM01V1___htot[MAX_MULT]; */
/*     float           MM01V1___hs[MAX_MULT]; */

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

/*     TBranch        *b_MM03Y1___hhits; */
/*     TBranch        *b_MM03Y1___hch; */
/*     TBranch        *b_MM03Y1___ht; */
/*     TBranch        *b_MM03Y1___htot; */
/*     TBranch        *b_MM03Y1___hs; */

/*     TBranch        *b_MM03X1___hhits; */
/*     TBranch        *b_MM03X1___hch; */
/*     TBranch        *b_MM03X1___ht; */
/*     TBranch        *b_MM03X1___htot; */
/*     TBranch        *b_MM03X1___hs; */

/*     TBranch        *b_MM03U1___hhits; */
/*     TBranch        *b_MM03U1___hch; */
/*     TBranch        *b_MM03U1___ht; */
/*     TBranch        *b_MM03U1___htot; */
/*     TBranch        *b_MM03U1___hs; */

/*     TBranch        *b_MM03V1___hhits; */
/*     TBranch        *b_MM03V1___hch; */
/*     TBranch        *b_MM03V1___ht; */
/*     TBranch        *b_MM03V1___htot; */
/*     TBranch        *b_MM03V1___hs; */

/*     TBranch        *b_MM01Y1___hhits; */
/*     TBranch        *b_MM01Y1___hch; */
/*     TBranch        *b_MM01Y1___ht; */
/*     TBranch        *b_MM01Y1___htot; */
/*     TBranch        *b_MM01Y1___hs; */

/*     TBranch        *b_MM01X1___hhits; */
/*     TBranch        *b_MM01X1___hch; */
/*     TBranch        *b_MM01X1___ht; */
/*     TBranch        *b_MM01X1___htot; */
/*     TBranch        *b_MM01X1___hs; */

/*     TBranch        *b_MM01U1___hhits; */
/*     TBranch        *b_MM01U1___hch; */
/*     TBranch        *b_MM01U1___ht; */
/*     TBranch        *b_MM01U1___htot; */
/*     TBranch        *b_MM01U1___hs; */

/*     TBranch        *b_MM01V1___hhits; */
/*     TBranch        *b_MM01V1___hch; */
/*     TBranch        *b_MM01V1___ht; */
/*     TBranch        *b_MM01V1___htot; */
/*     TBranch        *b_MM01V1___hs; */

   TBranch        *b_PB01U1___chits,*b_PB01U1___cch,*b_PB01U1___ct;   
   TBranch        *b_PB01X1___chits,*b_PB01X1___cch,*b_PB01X1___ct;   
   TBranch        *b_PB02V1___chits,*b_PB02V1___cch,*b_PB02V1___ct;   

   TBranch        *b_PB03U1___chits,*b_PB03U1___cch,*b_PB03U1___ct;   
   TBranch        *b_PB03X1___chits,*b_PB03X1___cch,*b_PB03X1___ct;   
   TBranch        *b_PB04V1___chits,*b_PB04V1___cch,*b_PB04V1___ct;   

   TBranch        *b_PB05U1___chits,*b_PB05U1___cch,*b_PB05U1___ct;   
   TBranch        *b_PB05X1___chits,*b_PB05X1___cch,*b_PB05X1___ct;   
   TBranch        *b_PB06V1___chits,*b_PB06V1___cch,*b_PB06V1___ct;   




   MMan(TFile *file, TTree *tree, bool fHitsInTree=true, int mult=0, float chi2=0);

   // will chain all chunks for this run 
   MMan(const char* dir,int run, const char* detdat, bool fHitsInTree=true, int mult=0, float chi2=0);

   ~MMan();

   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop(int nentries);
   Bool_t Notify();
   void   Show(Int_t entry = -1);

   void Allocate(int run);
   void TreatEvent();
   void Entry(Int_t i);
   void Test();

   void Efficiency(const char* chamber,const char* outfile,const char* opt); 
   void PlotTracking(int n);

   void EffVsChi2(const char* outfile, float chi2min, float chi2max, float ste);

   void TreeOn(bool dotree) {fDoTree = dotree;}
   
   // find time cuts for all planes (n sigmas)
   void TimeCuts(int n);

   ClassDef(MMan,0)
};



#endif









