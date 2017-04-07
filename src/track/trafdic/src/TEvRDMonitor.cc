// $Id: TEvRDMonitor.cc 14094 2015-11-06 15:28:48Z lsilva $

#include <iostream>
#include "TMath.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
//#define RDMon_AligTree
//#define RDMon_SM2Tree
#if defined RDMon_AligTree || defined RDMon_SM2Tree
#  include "TTree.h"
#endif
#include "CsInit.h"
//#define RDMon_GEMWrBack
#ifdef RDMon_GEMWrBack
#  include "CsOpt.h"
#endif
#include "CsEvent.h"
#include "CsHistograms.h"
#include "CsMCUtils.h"
#include "CsCluster.h"
#include "CsDigit.h"
#include "CsDetector.h"
#include "CsGEMDetector.h"
#include "CsPixelGEMDetector.h"
#include "CsRandom.h"
#include "DaqDataDecoding/DaqOption.h"
#include "Traffic.h"
#include "TEv.h"
#include "TSetup.h"
#include "TOpt.h"

using namespace CS;
using namespace std;

/*! 
  \brief Traffic method for monitoring RD reconstruction.
*/

/*
  Various residuals histogrammed depending upon "TOpt::Hist[1]">1 and
  "TOpt::Hist[16-18]"
 */

class THlx0 {       // Initial helix
public:
  THlx0();
  THlx0(const THlx &h): hlx(h) {};
  const THlx &hlx;
  bool operator<( const THlx0 &hlx0 ) const {  //!< "less than" operator
    if (hlx(0)<hlx0.hlx(0)) return true;
    else                    return false;
  }
};

void TEv::RDMonitor()
{
  if(TOpt::Hist[1] == 0) return;
  //if( IsMC() ) return;
  
  const TSetup& setup = TSetup::Ref();
  CsMCUtils utils;

  //----------------------------------------

  static int igr_mx;

  const double trackTMx[5] = // Max track "sigmaTime" allowed for sampling dets
    // 2.5 = 1 sigma of 12uM = 9ns / sqrt(12)
    { 5,2.5,2.5, 5,2.5};   // ...Non-drifts, DC, ST, DW, MB

  //    ******************** HISTOGRAMS ********************

  //     *************** MAIN HISTOGRAMS ***************
  static TH1D *h_trigger;
  // 4 copies: I(nclusive), S(emi)I(inclusive),
  //           C(alorimeter)(+, possibly, high-Q2), O(uter)
  static TH2D *h_nhits[4], *h_chi2[4], *h_prob[4], *h_mom[4];
  static TH1D *h_time[4], *h_tDt[4], *h_nBeams[4], *h_nBoTs[4];
  static vector<CsHist2D*> residus[2];
  static vector<CsHist2D*> reperps[2];
  static vector<CsHist1D*> guess, effic, rndom, bckgd;
  static int timingRequired = -1; // Flag used by the special histogramming section (residuals et sqq). !=0 means required. Default (=-1) is to be overwritten in the booking block.

  //    *************** AUXILLIARY HISTOGRAMS ***************
  // They are brought into play by a series of cpp macros...
  // ...Possibly nested.

  //#define RDMon_DISCARD_HALOS  // Discard halos, even when probing detectors downstream of muWall#2

  //#define RDMon_EFF_MAPs
#ifdef RDMon_EFF_MAPs
  static vector<CsHist2D*> guess2, effic2;
#endif

  //#define RDMon_DBL_LAYER
#ifdef RDMon_DBL_LAYER
  static map<int,int> dblhists;
  static vector<CsHist2D*> dblresol, dblreso2;
  //#  define RDMon_TOP_SPLIT_STRAWS 2  // Restrict to split straws (1:top, 2:all but top; 3...)
  static vector<CsHist1D*> dblguess, dbleffic, dblrndom, dblbckgd;
#endif

  //#define RDMon_ScifiSi_TIME	// Timing of scifis and Sis

  //#define RDMon_MMAmplitude	// MM amplitude per track (in view of dE/dx based PID)
#ifdef RDMon_MMAmplitude
  static TH2D *mMAvsP;
  static int mMAWord; static unsigned int mMAPats[2];
#endif

  //#define RDMon_DEBUG_MAs 1 // I.e. determine extent and position of dead zones. >1 for dumping events

  // ********** DIGIT DATA **********

  //#define RDMon_DIGIT_DATA 1	// =2: Allow also off-time tracks: the idea being that timing may not matter if one is interested in some particuliar digit data
#ifdef RDMon_DIGIT_DATA
  static map<int,int> mmhists, drifthists, gemhists;
  //#  define RDMon_DIGIT_MM
#  ifdef RDMon_DIGIT_MM
  static vector<CsHist1D*> mmltime;
  static vector<CsHist1D*> mmttime;
  static vector<CsHist1D*> mmctime;
  static vector<CsHist2D*> mmhtime;
  static vector<CsHist1D*> mmcToT;
  static vector<CsHist2D*> mmhtvT;   // ToT vs. t
#  endif
  //#  define RDMon_DIGIT_DC
#  ifdef RDMon_DIGIT_DC
  static vector<CsHist2D*> drifttime;
  static vector<CsHist2D*> driftres;
  static vector<CsHist1D*> driftges, drifteff, driftrdm, driftbck;
#  endif
#  define RDMon_DIGIT_GM    // GM and SI
#  define RDMon_nGMs_MX 74
#  ifdef RDMon_DIGIT_GM
#    define RDMon_GM_2D_CALIB
#    ifdef RDMon_GM_2D_CALIB
  static vector<CsHist2D*> gemampls[3];// HIT amplitude distributions: 3 samples
  static vector<CsHist2D*> gemnoise;   // ``Noise'' HIT amplitude distribution
#    else
  static vector<CsHist1D*> gemampls[3];// HIT amplitude distributions: 3 samples
  static vector<CsHist1D*> gemnoise;   // ``Noise'' HIT amplitude distribution
#    endif
  static vector<CsHist2D*> gemratio;   // HIT amplitude ratios: 3/2 vs. 2/1
  static vector<CsHist2D*> gemacorr;   // CLUSTER amplitude correlation
  static vector<CsHist1D*> gemadiff;   // CLUSTER amplitude diff
  static vector<CsHist1D*> gemctime;   // CLUSTER time
  // Calibration (intercept, slope, quad) of amplitude correlation X(Y)(V(U))...
  static vector<const float*> gemcorrs;
  // ...w/ default (0,1,0)...
  float gemdflt[3]; gemdflt[0] = 0; gemdflt[1] = 1; gemdflt[2] = 0;
  // ...Default is automatically assigned to all GEM if RDMon_GM_DEFAULT
  //#    define RDMon_GM_DEFAULT
#  endif
#endif

  //#define RDmon_SCIFI_TIME
#ifdef RDmon_SCIFI_TIME
  static CsHist1D *scifihists[5]; 
#endif

  // *************** ROOT Trees ***************

#ifdef RDMon_AligTree
  static TTree *AligTree;
  static float txAlig, tyAlig;  // TTrack angles
  static int nhAlig;        // # of all Alig hits
  static int nAlig[4];      // # of Alig hits in 4 zones (<SM1,<RICH,<SM2,>SM2)
  static unsigned int iAlig[2]; // Alig hits pattern
  static float uAlig[52];
  static vector<int> Alig_map;
  static int Alig_ipl0;
#endif
#ifdef RDMon_GEMWrBack
  static int GEMs_stream = -1;
#endif
#ifdef RDMon_SM2Tree
  static TTree *sm2_tree;
  static int   sm2_nhits;
#  define SM2_NPATS 5
  static unsigned int sm2_pats[SM2_NPATS];   // Valid for 2001 setup: 141 det's
  static float sm2_hlx[5];
  static float sm2_chi2;
  static unsigned int sm2_pguess, sm2_pfound;
  static float sm2_uguess[12], sm2_vguess[12], sm2_found[12];
#endif

  // INTERCEPTS of ttracks @ ORIGIN (0x3 and 0x6)
  // and @ beginning of Zone 1 (0x3) and 2 (0x6)
  static CsHist2D *intrcpts[6];
  static CsHist1D *vsslat;   // Type of TTrack according to VS, S or LAT
  //     ********** SINGLE OUT SPECIAL TYPES of TRIGGERS **********
  // - Special types are:
  //   - "stdTrigs" = standard triggers.
  //   - "specialTrigs", to be inclusively histo'd,
  //   - "pureTrigs", to be histo'd only when pure (typically their overlap w/
  //    the rest of the triggers is large).
  //  These 3 types allow to single out interesting events for histo'ing,
  //  whatever the data taking.
  // - Muon data taking: "specialTrigs" is assigned OT. "pureTrigs" is assigned
  //  pure CT (w/ allowed combination w/ highQ2T). And "stdTrigs" are IMLT.
  // - Hadron data taking: "specialTrigs" is assigned NIT, i.e. non interacting
  //  (interesting for it allows to monitor the reconstruction of paraxial
  //  tracks, much akin to that of the scattered hadron, but which is
  //  unfortunately not enabled in physics data taking). "pureTrigs" are
  //  beam-based triggers, BT|aBT (which, w/ their large overlap w/ the rest of
  //  the triggers, are interesting to consider only when pure). And "stdTrigs"
  //  are RPD-based triggers.
  // - DVCS running: "specialTrigs" is assigned "DVCS" in 2009, and "Tiger" in 
  //  2012. And "stdTrigs" are MLO[L]T. "pureTrigs" are BT|aBT.
  // - Primakoff: "stdTrigs" is assigned Primakoff1|2. "pureTrigs" are BT|aBT.
  //  "specialTrigs" is RPD.
  // - Low-t: "stdTrigs" is assigned lowt1|2. "pureTrigs" are BT|aBT.
  //  "specialTrigs" is RPD.
  // - DY: "stdTrigs" is assigned MTLT|OTLT|LAST2, "pureTrigs" is CT,
  //  "specialTrigs" is "MT|OT|LT".
  static unsigned int pureTrigs, highQ2Trig = 0, specialTrigs, stdTrigs = 0x7;

  static bool first = true; if (first) {
    first = false; // ******************** INITIALISATIONS ********************

    bool hadronRun = CsInit::Instance()->IsAHadronJob();
    bool primakoffRun, lowtRun = false, dyRun; int  dvcsRun;

    const unsigned int caloTrig = 0x10;
    unsigned int caloBasedTrigs = caloTrig, beamTrigs = 0;
    const CS::DaqOption &opts = CsInit::Instance()->getDaqEventsManager().GetDaqOptions();
    int bit = -1;// Bit 0x200 doesn't necessarily mean high-Q2 trigger: check it
    try { bit = opts.GetTriggerBit("LargeQ2Trigger"); } catch ( ... ) { }
    if (bit>=0) { highQ2Trig = 1<<bit; caloBasedTrigs |= highQ2Trig; }
    bit = -1;
    try { bit = opts.GetTriggerBit("BeamTrigger"); } catch ( ... ) { }
    if (bit>=0) beamTrigs |= 1<<bit;
    bit = -1;
    try { bit = opts.GetTriggerBit("altBeamTrigger"); } catch ( ... ) { }
    if (bit>=0) beamTrigs |= 1<<bit;
    pureTrigs =  hadronRun ? beamTrigs : caloBasedTrigs;
    unsigned int niTrig =// As of 2008/12, NIT's not described in the mapping...
      0x10; // ...=> Built-in (although NIT is not present throughout)
    specialTrigs = hadronRun ? niTrig : 0x8;
    // ***** DVCS data taking
    bit = -1;
    try { bit = opts.GetTriggerBit("DVCS"); } catch ( ... ) { }
    if (bit>=0) {
      specialTrigs = 1<<bit; pureTrigs = beamTrigs; stdTrigs = 0xe;
      dvcsRun = 2009;
    }
    else {
      dvcsRun = 0;
    }
    try { bit = opts.GetTriggerBit("Tiger"); } catch ( ... ) { }
    if (bit>=0) {
      specialTrigs = 1<<bit; pureTrigs = beamTrigs; stdTrigs = 0x20e;
      dvcsRun = 2012;
    }
    // ***** DY data taking
    bit = -1;
    try { bit = opts.GetTriggerBit("LAST2"); } catch ( ... ) { }
    if (bit>=0) {
      specialTrigs = 0x20a; pureTrigs = 0x10; stdTrigs = 0x105;
      dyRun = true;
    }
    else {
      dyRun = false;
    }
    // ***** Primakoff data taking
    bit = -1;
    try { bit = opts.GetTriggerBit("Primakoff1"); } catch ( ... ) { }
    if (bit<0) {
      try { bit = opts.GetTriggerBit("Prim1"); } catch ( ... ) { }
    }
    if (bit>=0) {
      stdTrigs = 1<<bit; primakoffRun = true;
    }
    else primakoffRun = false;
    bit = -1;
    try { bit = opts.GetTriggerBit("Primakoff2"); } catch ( ... ) { }
    if (bit<0) {
      try { bit = opts.GetTriggerBit("Prim2"); } catch ( ... ) { }
    }
    if (bit>=0) {
      if (primakoffRun) stdTrigs |= 1<<bit;
      else              stdTrigs = 1<<bit;
      primakoffRun = true;
    }
    if (primakoffRun) {
      pureTrigs = beamTrigs; specialTrigs = 0x1;
    }
    else {
      // ***** Low-t data taking
      // (Note: This is conditioned by !=Primakoff, because lowt triggers are
      // still described in the mapping (DAQ.xml) for the 2009 Primakoff data
      // taking, although they, probably, were not active.)
      bit = -1;
      try { bit = opts.GetTriggerBit("lowt1"); } catch ( ... ) { }
      if (bit>=0) {
	stdTrigs = 1<<bit; lowtRun = true;
      }
      else lowtRun = false;
      bit = -1;
      try { bit = opts.GetTriggerBit("lowt2"); } catch ( ... ) { }
      if (bit>=0) {
	if (lowtRun) stdTrigs |= 1<<bit;
	else         stdTrigs = 1<<bit;
	lowtRun = true;
      }
      if (lowtRun) {
	pureTrigs = beamTrigs; specialTrigs = 0x1;
      }
    }

    //           *************** HISTO BOOKING ***************
    CsHistograms::SetCurrentPath("/Traffic/RDmonitor");

    //         ********** TRACKS (as function of Trigger Type) **********
    h_trigger = new TH1D("h_trigger","Trigger",4096,-.5,4095.5);
    char hN[]  = "si_nbeams"; char hT[] = "#chi^{2} prob - PrimakoffT   ";
    const char *tabs[4] = {"i","si","c","o"};
    const char *tags[7][4]= { {""," - ILMT",      " - CT",   " - OT"},
			      {""," - RPDT",      " - (a)BT"," - NIT"},
			      {""," - LMOT",      " - (a)BT"," - DVCST"},
			      {""," - PrimakoffT"," - (a)BT"," - DT0"},
			      {""," - lowtT",     " - (a)BT"," - DT0"},
			      {""," - LMOLT",     " - (a)BT"," - TigerT"},
			      {""," - 2muT",      " - CT",   " - 1muT"} };
    int i = 0; if (hadronRun) i = 1;
    if (dvcsRun==2009) i = 2; if (dvcsRun==2012) i = 5; if (dyRun) i = 6;
    if (primakoffRun) i = 3; if (lowtRun) i = 4;

    for (int j = 0; j<4; j++) {
      sprintf(hN,"%s_nhits",tabs[j]);  sprintf(hT,"#Hits%s",        tags[i][j]);
      h_nhits[j] = new TH2D(hN,hT,120,  0 ,120,16,-.5,15.5 );
      sprintf(hN,"%s_chi2",tabs[j]);   sprintf(hT,"#chi^{2}/ndf%s", tags[i][j]);
      h_chi2[j]  = new TH2D(hN,hT, 50,  0 , 30,16,-.5,15.5);
      sprintf(hN,"%s_prob",tabs[j]);   sprintf(hT,"#chi^{2} prob%s",tags[i][j]);
      h_prob[j]  = new TH2D(hN,hT,100,  0 ,  1,16,-.5,15.5);
      sprintf(hN,"%s_mom", tabs[j]);   sprintf(hT,"P (GeV)%s",      tags[i][j]);
      double pMn, pMx; if (TOpt::iCut[15]>0) { pMn = -100; pMx = 250; }
      else                                   { pMn = -250; pMx = 100; }
      h_mom[j]   = new TH2D(hN,hT,175,pMn,pMx, 16,-.5,15.5);
      sprintf(hN,"%s_time",tabs[j]);   sprintf(hT,"Time%s",         tags[i][j]);
      h_time[j]  = new TH1D(hN,hT,100,-25,25);
      sprintf(hN,"%s_tDt", tabs[j]);   sprintf(hT,"t/#deltat%s",    tags[i][j]);
      h_tDt[j]   = new TH1D(hN,hT,100,-50,50);
      sprintf(hN,"%s_nBeams",tabs[j]); sprintf(hT,"#scifi/Si's%s",  tags[i][j]);
      h_nBeams[j]= new TH1D(hN,hT,16,-0.5,15.5);
      sprintf(hN,"%s_nBoTs",tabs[j]);  sprintf(hT,"#BoT's%s",       tags[i][j]);
      h_nBoTs[j] = new TH1D(hN,hT,16,-0.5,15.5);
    }
#ifdef RDMon_MMAmplitude
    mMAvsP = new TH2D("mMAvsP","MM Amplitude vs. P",100,0,40,50,0,2.5);
    int idet; // Fill pattern of MM planes
    for (idet = setup.vIplFirst()[0], mMAWord = -1, mMAPats[0]=mMAPats[1] = 0;
	 idet<=setup.vIplLast()[0]; idet++) {
      const TDetect &d = setup.vDetect()[idet]; if (d.IType!=27) continue;
      if (mMAWord==-1) mMAWord = idet/32;
      mMAPats[idet/32-mMAWord] |= 1<<idet%32;
    }
#endif

    if (TOpt::Hist[1]&0x10) {

      // ****************************** RESIDUALS ******************************


#ifdef RDMon_AligTree
      AligTree = new TTree("T","tracks");
      AligTree->Branch("txAlig",&txAlig,"txAlig/F",32000);
      AligTree->Branch("tyAlig",&tyAlig,"tyAlig/F",32000);
      AligTree->Branch("nhAlig",&nhAlig,"nhAlig/I",32000);
      AligTree->Branch("nAlig" ,nAlig,  "nAlig[4]/I",32000);
      AligTree->Branch("uAlig" ,uAlig,  "uAlig[nhAlig]/F",32000);
      AligTree->Branch("iAlig" ,iAlig,  "iAlig[2]/i",32000); // i: unsigned
      nhAlig=nAlig[0]=nAlig[1]=nAlig[2]=nAlig[3] = 0; iAlig[0]=iAlig[1] = 0;
      Alig_ipl0 = -1;
      txAlig=tyAlig = 0;
#endif
#ifdef RDMon_GEMWrBack
      string filename;
      if (CsOpt::Instance()->getOptRec("RDMon_GEMWrBack","file",filename))
	GEMs_stream = CsEvent::Instance()->openRawEventOutputStream(filename);
#endif
#ifdef RDMon_SM2Tree
      sm2_tree = new TTree("S","sm2_tracks");
      sm2_tree->Branch("nhits"    ,&sm2_nhits, "nhits/I"       ,32000);
      sm2_tree->Branch("dets_pats",&sm2_pats,  "dets_pats[5]/i",32000);
      sm2_tree->Branch("helix"    ,&sm2_hlx,   "helix[5]/F"    ,32000);
      sm2_tree->Branch("chi2"     ,&sm2_chi2,  "chi2/F"        ,32000);
      sm2_tree->Branch("guess_pat",&sm2_pguess,"guess_pat/i"   ,32000);
      sm2_tree->Branch("found_pat",&sm2_pfound,"found_pat/i"   ,32000);
      sm2_tree->Branch("uguess"   ,&sm2_uguess,"uguess[12]/F"  ,32000);
      sm2_tree->Branch("vguess"   ,&sm2_vguess,"vguess[12]/F"  ,32000);
      sm2_tree->Branch("found"    ,&sm2_found, "found[12]/F"   ,32000);
#endif
#ifdef RDMon_DBL_LAYER
      int idbl= 0;
#endif
#ifdef RDMon_DIGIT_DATA
      int imm = 0, idrift = 0, igem = 0;
#endif
      int ipl;
      // Arrays of char's for histo titles: to be edited w/ sprintf
      // ...CAUTION! Check format strings fit in allocated length
      char htitle[] = "DC01V2__ Residuals Zone 0 vs. V+#pi/2";
      igr_mx = (int)setup.vIplFirst().size();
      for (int igr = 0; igr<igr_mx; igr++) {
	if ((1<<igr&(TOpt::Hist[17]|   // Planes are to be monit'd if inactive
		     TOpt::Hist[16]))  // Planes are to be monitored whatsoever
	    ==0) continue;

	// *************** BOOK per plane HISTOs ***************

	for (ipl = setup.vIplFirst()[igr]; ipl<=setup.vIplLast()[igr]; ipl++) {

	  // ********** LOOP OVER PLANES IN ZONE "igr" **********

	  const TPlane &p = setup.vPlane(ipl);
	  const TDetect &d = setup.vDetect(p.IDetRef);

#ifdef RDMon_DIGIT_GM
	  // ***** SPECIAL CASE of GEM/Si AMPLITUDES... *****
	  if (d.IType!=21 && d.IType!=26 && d.IType!=28) continue;
#else
	  if (p.IFlag && (1<<igr&TOpt::Hist[16])==0) continue;
#endif

#ifdef RDMon_AligTree
	  // Build map: ipl (-ipl0) -> igem
	  if (igr<3) {
	    if (Alig_ipl0<0) Alig_ipl0 = ipl;
	    if (d.IType==26 || d.IType==27) {  // GMs or MMs
	      if (nhAlig==52) {
		cout << "# GMs + # MMs > 52!\n"; exit(1);
	      }
	      // List of wirD's (to be stored in 1st TTree event)
	      uAlig[nhAlig] = d.PtrDet()->getWirD();
	      iAlig[nhAlig/32] |= 1<<nhAlig%32;  // Dets pattern
	      Alig_map.push_back(nhAlig++);      // Build map of Alig dets
	      if      (d.X(0)< 350) nAlig[0]++;  // Count dets in 4 zones
	      else if (d.X(0)< 750) nAlig[1]++;
	      else if (d.X(0)<1750) nAlig[2]++;
	      else                  nAlig[3]++;
	    }
	    else
	      Alig_map.push_back(-1);
	  }
#endif

	  // ***** HIGH RESOLUTION: residus[0], reperps[0] *****
	  double delta;
	  if (d.PtrDet()->hasDrift()) delta = d.Pitch/2;
	  else if (d.Name.find("SI")==0)  // Since det is much more precise...
	    /* ...than tracking */    delta = d.Pitch*10;
	  else                        delta = d.Pitch*3;
	  //#define VERY_LOW_RESOLUTION 2
#ifdef VERY_LOW_RESOLUTION
	  delta *= 4*VERY_LOW_RESOLUTION;
#endif

	  if ((d.IType==39 || d.IType==40) && timingRequired==-1)
	    timingRequired = 0; // Do not require timing when histogramming MAs and MAs only. Because these detectors are expected to be fired mainly by pure calo trigger, for which timing is problematic. And since their timing does not matter so much...
	  else timingRequired = 1;

	  char coord = d.Name.c_str()[4];
	  int nbins = d.Nwires<50 ? d.Nwires : 25;
	  //if (strncmp(d.Name.c_str(),"ST",2)==0) nbins = d.Nwires;

	  sprintf(htitle,"%s Residuals Zone 0 vs. %c",d.Name.c_str(),coord);
	  double C0 = d.X(1)*d.Ca+d.X(2)*d.Sa;
	  double CSiz = fabs(d.Ca)>.02 ? d.Siz(1) : d.Siz(2);
	  residus[0].push_back
	    (new CsHist2D(d.Name+"_res0",htitle, 
			  nbins,C0-CSiz/2,C0+CSiz/2,150,-delta,delta));
	  sprintf(htitle,"%s Residuals Zone 0 vs. %c+#pi/2",d.Name.c_str(),coord);
	  double P0 = d.X(2)*d.Ca-d.X(1)*d.Sa;
	  double PSiz = fabs(d.Ca)>.02 ? d.Siz(2) : d.Siz(1);
	  reperps[0].push_back
	    (new CsHist2D(d.Name+"_rep0",htitle, 
			  nbins,P0-PSiz/2,P0+PSiz/2,150,-delta,delta));
	  // ***** LOW RESOLUTION: residus[1], reperps[1] *****
	  if (d.PtrDet()->hasDrift()) {
	    //#define RDMon_DW_IMPRECISE
#ifdef RDMon_DW_IMPRECISE
	    // Imprecise tracking (To be enabled in 2002, where DWs are isolated
	    if (d.Name.find("DW")==0) delta = d.Pitch*1.5;
#endif
	    delta = d.Pitch;
	  }
	  else if (d.Name.find("SI")==0)  // Since det is much more precise...
	    /* ...than tracking */    delta = d.Pitch*50;
	  else if (d.Name.find("H")==0)   // Since det is much less precise...
	    /* ...than tracking */    delta = d.Pitch*5;
	  else if (d.Name.find("PA")==0)  // Since det is close to SAT
	    /* ...              */    delta = d.Pitch*5;
	  else if (d.Name.find("MB")==0 && d.Name.c_str()[7]=='c')  // 
	    /* ...              */    delta = d.Pitch*40;
	  else                        delta = d.Pitch*10;
#ifdef VERY_LOW_RESOLUTION
	  delta *= 4*VERY_LOW_RESOLUTION;
#endif
	  sprintf(htitle,"%s Residuals Zone 1 vs. %c",d.Name.c_str(),coord);
	  residus[1].push_back
	    (new CsHist2D(d.Name+"_res1",htitle, 
			  nbins,C0-CSiz/2,C0+CSiz/2,150,-delta,delta));
	  sprintf(htitle,"%s Residuals Zone 1 vs. %c+#pi/2",d.Name.c_str(),coord);
	  reperps[1].push_back
	    (new CsHist2D(d.Name+"_rep1",htitle, 
			  nbins,P0-PSiz/2,P0+PSiz/2,150,-delta,delta));

	  // ********** EFFICIENCY: "guess' and "effic"... **********
	  sprintf(htitle,"%s Guess vs. %c",d.Name.c_str(),coord);
	  guess.push_back
	    (new CsHist1D(d.Name+"_guess",htitle, 
			  d.Nwires,C0-CSiz/2,C0+CSiz/2));
	  sprintf(htitle,"%s Effic vs. %c",d.Name.c_str(),coord);
	  effic.push_back
	    (new CsHist1D(d.Name+"_effic",htitle, 
			  d.Nwires,C0-CSiz/2,C0+CSiz/2));
	  sprintf(htitle,"%s Random vs. %c",d.Name.c_str(),coord);
	  rndom.push_back
	    (new CsHist1D(d.Name+"_rndom",htitle, 
			  d.Nwires,C0-CSiz/2,C0+CSiz/2));
	  sprintf(htitle,"%s Background vs. %c",d.Name.c_str(),coord);
	  bckgd.push_back
	    (new CsHist1D(d.Name+"_bckgd",htitle, 
			  d.Nwires,C0-CSiz/2,C0+CSiz/2));
#ifdef RDMon_EFF_MAPs
	  // ... and 2D maps
#  define RDMon_MAP_SPECIAL 1
#  ifdef RDMon_MAP_SPECIAL
	  if (d.Name.find("H")==0 ||
	      d.Name.find("MB")==0 ||
	      d.Name.find("MA")==0 ||
	      d.Name.find("DW")==0) {
	    int gran = RDMon_MAP_SPECIAL;
	    sprintf(htitle,"%s Guess 2D",d.Name.c_str());
	    guess2.push_back
	      (new CsHist2D(d.Name+"_guess2",htitle, 
			    d.Nwires*gran,C0-CSiz/2,C0+CSiz/2,
			    d.Nwires*gran,P0-PSiz/2,P0+PSiz/2));
	    sprintf(htitle,"%s Effic 2D",d.Name.c_str());
	    effic2.push_back
	      (new CsHist2D(d.Name+"_effic2",htitle, 
			    d.Nwires*gran,C0-CSiz/2,C0+CSiz/2,
			    d.Nwires*gran,P0-PSiz/2,P0+PSiz/2));
	  }
	  else
#  endif
	    {
	      int Nbins;
	      if (d.PtrDet()->hasDrift()) Nbins = d.Nwires/8;
	      else                        Nbins = d.Nwires/16;
	      sprintf(htitle,"%s Guess 2D",d.Name.c_str());
	      guess2.push_back
		(new CsHist2D(d.Name+"_guess2",htitle, 
			      Nbins,C0-CSiz/2,C0+CSiz/2,
			      Nbins,P0-PSiz/2,P0+PSiz/2));
	      sprintf(htitle,"%s Effic 2D",d.Name.c_str());
	      effic2.push_back
		(new CsHist2D(d.Name+"_effic2",htitle, 
			      Nbins,C0-CSiz/2,C0+CSiz/2,
			      Nbins,P0-PSiz/2,P0+PSiz/2));
	    }
#endif
#ifdef RDMon_DBL_LAYER
	  if (strncmp(d.Name.c_str(),"DC",2)==0 ||
	      strncmp(d.Name.c_str(),"DW",2)==0 ||
	      strncmp(d.Name.c_str(),"ST",2)==0
	      /* && d.Name[7]=='b' */   ||              // "b" slice
	      strncmp(d.Name.c_str(),"MB",2)==0 
	      /* && (d.Name[7]=='b' || d.Name.[7]=='r') */) {
	    dblhists[d.IDet] = idbl++;
	    static float dblxprv = 0;
	    double dx_mx = 1.5;  // Max distance between 2 layers of DL
	    if      (d.Name[0]=='M') dx_mx = 3;
	    else if (d.Name[1]=='W') dx_mx = 2.5;
	    if (fabs(d.X(0)-dblxprv)<dx_mx) {
	      // Current det close (<1cm) to previous: planes of same detector
	      // => book double layer histograms
	      char det_name[] = "ST03X1bud";
	      char station = d.Name[3], coord = d.Name[4];
	      if (strncmp(d.Name.c_str(),"DC",2)==0) {
		sprintf(det_name,"DC0%c%c12",station,coord);
	      }
	      else if (strncmp(d.Name.c_str(),"ST",2)==0) {
		char number = d.Name[5], slice = d.Name[7];
		sprintf(det_name,"ST0%c%c%c%cud",station,coord,number,slice);
	      }
	      else if (strncmp(d.Name.c_str(),"MB",2)==0) {
		char number = d.Name[5], slice = d.Name[7];
		sprintf(det_name,"MB0%c%c%c%cud",station,coord,number,slice);
	      }
	      else if (strncmp(d.Name.c_str(),"DW",2)==0) {
		sprintf(det_name,"DW0%c%c12",station,coord);
	      }
	      sprintf(htitle,"%s RDM DL Resol",det_name);
	      dblresol.push_back
		(new CsHist2D((string)det_name+"_dlres",
			      htitle,150,-d.Pitch/2,d.Pitch/2,
			      3,C0-CSiz/2,C0+CSiz/2));
	      dblreso2.push_back
		(new CsHist2D((string)det_name+"_dlre2",
			      htitle,150,-d.Pitch/2,d.Pitch/2,
			      3,C0-CSiz/2,C0+CSiz/2));
	      sprintf(htitle,"%s RDM DL Guess",det_name);
	      dblguess.push_back
		(new CsHist1D((string)det_name+"_dlges",
			      htitle,d.Nwires,C0-CSiz/2,C0+CSiz/2));
	      sprintf(htitle,"%s RDM DL Effic",det_name);
	      dbleffic.push_back
		(new CsHist1D((string)det_name+"_dleff",
			      htitle,d.Nwires,C0-CSiz/2,C0+CSiz/2));
	      sprintf(htitle,"%s RDM DL Rndom",det_name);
	      dblrndom.push_back
		(new CsHist1D((string)det_name+"_dlrdm",
			      htitle,d.Nwires,C0-CSiz/2,C0+CSiz/2));
	      sprintf(htitle,"%s RDM DL Bckgd",det_name);
	      dblbckgd.push_back
		(new CsHist1D((string)det_name+"_dlbck",
			      htitle,d.Nwires,C0-CSiz/2,C0+CSiz/2));
	    }
	    dblxprv = d.X(0);
	  }
#endif
#ifdef RDMon_DIGIT_DATA
	  // ********** DIGIT DATA **********
#  ifdef RDMon_DIGIT_MM
	  if (strncmp(d.Name.c_str(),"MM",2)==0) {
	    mmhists[d.IDet] = imm++;
	    sprintf(htitle,"%s RDM Lead T",d.Name.c_str());
	    mmltime.push_back
	      (new CsHist1D(d.Name+"_lt",htitle,150,-325,390));
	    sprintf(htitle,"%s RDM Trail T",d.Name.c_str());
	    mmttime.push_back
	      (new CsHist1D(d.Name+"_tt",htitle,150,-325,390));
	    sprintf(htitle,"%s RDM Cluster T",d.Name.c_str());
	    mmctime.push_back
	      (new CsHist1D(d.Name+"_ct",htitle,150,-100,100));
	    sprintf(htitle,"%s RDM Hit T",d.Name.c_str());
	    mmhtime.push_back
	      (new CsHist2D(d.Name+"_ht",htitle,d.Nwires,-.5,d.Nwires-.5
			    ,150,-100,100));
	    sprintf(htitle,"%s RDM Cluster ToT",d.Name.c_str());
	    mmcToT.push_back
	      (new CsHist1D(d.Name+"_ctot",htitle,150,0,400));
	    sprintf(htitle,"%s RDM cToT vs. ct",d.Name.c_str());
	    mmhtvT.push_back
	      (new CsHist2D(d.Name+"_htvT",htitle,100,0,400,100,-100,100));
	  }
#  endif
#  ifdef RDMon_DIGIT_DC
	  if (strncmp(d.Name.c_str(),"ST",2)==0 ||
	      strncmp(d.Name.c_str(),"DC",2)==0 ||
	      strncmp(d.Name.c_str(),"DW",2)==0 ||
	      strncmp(d.Name.c_str(),"MB",2)==0) {
	    drifthists[d.IDet] = idrift++;
	    sprintf(htitle,"%s RDM dTime",d.Name.c_str());
	    drifttime.push_back
	      (new CsHist2D(d.Name+"_dt",htitle,32,0,d.PtrDet()->getTGate(),
			    50,-d.Pitch/2,d.Pitch/2));
	    sprintf(htitle,"%s RDM dTime",d.Name.c_str());
	    driftres.push_back
	      (new CsHist2D(d.Name+"_res",htitle,20,-d.Pitch/2,d.Pitch/2,
			    50,-d.Pitch/4,d.Pitch/4));
	    sprintf(htitle,"%s RDM GUESS",d.Name.c_str());
	    driftges.push_back
	      (new CsHist1D(d.Name+"_ges",htitle,20,-d.Pitch/2,d.Pitch/2));
	    sprintf(htitle,"%s RDM EFFIC",d.Name.c_str());
	    drifteff.push_back
	      (new CsHist1D(d.Name+"_eff",htitle,20,-d.Pitch/2,d.Pitch/2));
	    sprintf(htitle,"%s RDM RNDOM",d.Name.c_str());
	    driftrdm.push_back
	      (new CsHist1D(d.Name+"_rdm",htitle,20,-d.Pitch/2,d.Pitch/2));
	    sprintf(htitle,"%s RDM BCKGD",d.Name.c_str());
	    driftbck.push_back
	      (new CsHist1D(d.Name+"_bck",htitle,20,-d.Pitch/2,d.Pitch/2));
	  }
#  endif
#  ifdef RDMon_DIGIT_GM
	  if (strncmp(d.Name.c_str(),"GM",2)==0 ||
	      strncmp(d.Name.c_str(),"GP",2)==0 && d.Name[4]!='P' ||
	      strncmp(d.Name.c_str(),"SI",2)==0) {
	    if (igem>=RDMon_nGMs_MX)
	      CsErrLog::mes(elFatal,"RDMon_DIGIT_GM: Too many GEMs!");
	    int NWir = d.PtrDet()->getNWir();
	    gemhists[d.IDet] = igem++;
	    // HIT amplitude distributions (there are 3)
	    int iampl; char i[] = "0";
	    for (iampl = 0, i[0] = '0'; iampl<3; iampl++, i[0]++) {
	      sprintf(htitle,"%s RDM hAmpl%d",d.Name.c_str(),iampl);
#    ifdef RDMon_GM_2D_CALIB
	      gemampls[iampl].push_back
		(new CsHist2D(d.Name+"_ampl"+(string)i,htitle,
			      NWir,-.5,NWir-.5,64,0,1024));
#    else
	      gemampls[iampl].push_back
		(new CsHist1D(d.Name+"_ampl"+(string)i,htitle,128,0,1024));
#    endif
	    }
	    // HIT amplitude difference
	    sprintf(htitle,"%s RDM 2/0 vs. 1/0",d.Name.c_str());
	    gemratio.push_back
	      (new CsHist2D(d.Name+"_ratio",htitle,64,0,2,64,0,2));
	    // Attempt at catching noise
	    sprintf(htitle,"%s RDM hNoise",d.Name.c_str());
#    ifdef RDMon_GM_2D_CALIB
	    gemnoise.push_back
	      (new CsHist2D(d.Name+"_noise",htitle,
			    NWir,-.5,NWir-.5,32,0,512));
#    else
	    gemnoise.push_back
	      (new CsHist1D(d.Name+"_noise",htitle,64,0,512));
#    endif
	    static float gemxprv = 0;
	    if (fabs(d.X(0)-gemxprv)<.1) {
	      // Current det close (<1mm) to previous: planes of same detector
	      // => book CLUSTER amplitude correlations
	      char det_name[] = "GM01XY";
	      for (int c = 0; c<5; c++) det_name[c] = d.Name[c];
	      char coord = d.Name[4], counterpart;
	      if (coord=='X' || coord=='U') counterpart = coord+1;
	      else                          counterpart = coord-1;
	      det_name[5] = counterpart;
	      sprintf(htitle,"%s RDM a%c vs. a%c",det_name,counterpart,coord);
	      gemacorr.push_back
		(new CsHist2D((string)det_name+"_corr",
			      htitle,128,0,2048,128,0,2048));
	      sprintf(htitle,"%s RDM a%c - a%c",det_name,counterpart,coord);
	      gemadiff.push_back
		(new CsHist1D((string)det_name+"_diff",
			      htitle,250,-500,500));
#    ifdef RDMon_GM_DEFAULT  // Upon option, use default correlation coeff's...
	      // ...It's useful for determining new calibration coeff's.
	      gemcorrs.push_back(gemdflt);
#    else
	      if (d.IType==26 || d.IType==28) {
		CsGEMDetector *gd = dynamic_cast<CsGEMDetector*>(d.PtrDet());
		if (gd) 
		  gemcorrs.push_back(gd->getAmpCorr());
		else {
		  CsPixelGEMDetector *gpd =
		    dynamic_cast<CsPixelGEMDetector*>(d.PtrDet());
		  if (gpd) gemcorrs.push_back(gpd->getAmpCorr());
		  else
		    CsErrLog::msg(elFatal,__FILE__,__LINE__,
  "TDetect \"%s\" of type %d cannot be cast into either CsGEM or CsPixelGEM",
				  d.Name.c_str(),d.IType);
		}
	      }
	      else gemcorrs.push_back(gemdflt);
#    endif
	    }
	    sprintf(htitle,"%s RDM cTime",d.Name.c_str());
	    if (strncmp(d.Name.c_str(),"GM",2)==0)
	      gemctime.push_back
		(new CsHist1D(d.Name+"_ct",htitle,100,-50,50));
	    else
	      gemctime.push_back
		(new CsHist1D(d.Name+"_ct",htitle,100,-10,10));
	    gemxprv = d.X(0);
	  }
#  endif
#endif
	}   // End loop on planes

	// *************** BOOK per zone HISTOs ***************

#ifdef RDmon_SCIFI_TIME
	if (igr>=5) {
	  CsErrLog::mes(elFatal,"RDMon_SCIFI_TIME: Too many zones!");
	}
	char hname[] = "SciFiT_Z0";
	sprintf(hname,"SciFiT_Z%d",igr);
	sprintf(htitle,"Scifi Time Z%d",igr);
	scifihists[igr] = new CsHist1D(hname,htitle,100,-5,5);
#endif

      }   // End loop on group
#ifdef RDMon_AligTree
      AligTree->Fill();   // ***** Alig Tree: List of wirD's in 1st event *****
#endif

      if (TOpt::Hist[16]&0x20) {

	// ********** INTERCEPTS of TTracks @ ORIGINS **********
	// (booked with, dummy, zone=5)

	char hname[] = "YZ_X00_0x3"; char htitle[] = "YZ @ X==999. cm (0x3)";
	int iintrcpt, igr, igr_t;
	for (igr=iintrcpt = 0; igr<3 && igr<igr_mx; igr++) {
	  double X0;
	  if (igr==0) X0 = setup.TargetCenter[0];
	  else        X0 = setup.iPlane2Detect(setup.vIplFirst()[igr]).X(0);
	  if (igr==0 || igr==1) {
	    sprintf(hname,"YZ_X0%d_0x3",igr);
	    sprintf(htitle,"YZ @ X==%.0f cm (0x3)",X0);
	    if (igr==0) intrcpts[iintrcpt++] =
			  new CsHist2D(hname,htitle,50,-40,40,50,-40,40);
	    else        intrcpts[iintrcpt++] =
			  new CsHist2D(hname,htitle,100,-100,100,100,-100,100);

	  }
	  if (igr==0 || igr==2) {
	    sprintf(hname,"YZ_X0%d_0x6",igr);
	    sprintf(htitle,"YZ @ X==%.0f cm (0x6)",X0);
	    if (igr==0) intrcpts[iintrcpt++] =
			  new CsHist2D(hname,htitle,50,-40,40,50,-40,40);
	    else        intrcpts[iintrcpt++] =
			  new CsHist2D(hname,htitle,50,-80,80,50,-80,80);
	  }
	}
	intrcpts[iintrcpt++] =
	  new CsHist2D("YZ_MF2_gt","YZ @ X=38 m xx0>80",
		       50,-50,135,50,-50,50);
	intrcpts[iintrcpt++] =
	  new CsHist2D("YZ_MF2_le","YZ @ X=38 m xx0<=80",
		       50,-50,135,50,-50,50);
      }
      // Type of TTracks according to VS,S,LAT
      vsslat = new CsHist1D("vsslat","VS,S,LAT * 0x8*Zone",511,.5,511.5);

    }

    CsHistograms::SetCurrentPath("/"); // ********** END of BOOKING **********
  }  // *************** END of INITIALISATION BLOCK ***************

  //----------------------------------------
  
  if (CsEvent::Instance()->getTimeInSpill()<
      CsInit::Instance()->getMinTimeInSpill()) return;

  const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
  unsigned int evtTrig = trig_mask&allTrigs; // Cut away trailing bits

  {
    //#define DEBUG_TRACKING
#ifdef DEBUG_TRACKING
    static int ievt = 0, nsp1t = 0, nsp2t = 0, nmust = 0;
    ievt++; int nsp1 = 0, nsp2 = 0, nmus = 0;
#endif
    h_trigger->Fill(evtTrig);
    list<TTrack>::iterator it; int nBeams, nBoTs;
    for (it = listTrack.begin(), nBeams=nBoTs = 0; it != listTrack.end(); it++) {
      //   *************** HISTOs per TRACK: Loop over TTracks ***************
      TTrack &t = *it; THlx &H = t.Hfirst;
      if (t.Type&0x10) { nBeams++; if (t.Type&0x1) nBoTs++; }
      if (!H.with_mom() && (t.Type&0x10)==0) continue;
      // Fill histo's as a function of track's type (0 instead of 16 for Zone 5)
      for (int isic = 0; isic<4; isic++) {

	// Loop on Inclusive/Semi-Inclusive/Calorimeter
	if (isic==1 && !(evtTrig&stdTrigs) ||// Single out ``standard'' triggers
	    isic==2 &&
	    (evtTrig&pureTrigs)!=evtTrig ||  // Pure CT (+? highQ2) or beam
	    isic==3 && !(evtTrig&specialTrigs)) // OT or NIT or DVCST
	  continue;

	h_nhits[isic]->Fill(t.NHits+0.5,            t.Type%16);
	if (H.with_mom()) {
	  h_chi2[isic]->Fill(t.Chi2tot/(t.NDFs-5.), t.Type%16);
	  h_prob[isic]
	    ->Fill(TMath::Prob(t.Chi2tot,t.NDFs-5), t.Type%16);
	}
	else {
	  h_chi2[isic]->Fill(t.Chi2tot/(t.NDFs-4.), t.Type%16);
	  h_prob[isic]
	    ->Fill(TMath::Prob(t.Chi2tot,t.NDFs-4), t.Type%16);
	}
#ifdef DEBUG_TRACKING
	if ((t.Type&0x3)==0x3) nsp1++;
	if ((t.Type&0x6)==0x6) nsp2++;
	if ((t.Type&0xe)==0xe) nmus++;
#endif
	h_mom[isic]->Fill(1/H(5)                  , t.Type%16);
	h_time[isic]->Fill(t.MeanTime);
	h_tDt[isic]->Fill(t.MeanTime/t.SigmaTime);
      }
#ifdef RDMon_MMAmplitude
      if (fabs(H(5))>1/2.5) { // Fill MM amplitude only for low P tracks (where dE/dx based PID is expected to be possibly feasible)
	list<int>::const_iterator ih, ip; int nMMs; double amplitude;
	for (ih = t.lHPat().begin(), ip = t.lPRef().begin(), nMMs = 0,
	       amplitude = 0; ih!=t.lHPat().end(); ih++, ip++) {
	  if (*ih<0) continue;
	  int ipl = *ip, w = ipl/32-mMAWord; if (w<0) continue; if (1<w) break;
	  if (!(1<<ipl%32&mMAPats[w])) continue;
	  THit &h = vecHit[*ih];
	  double a; h.ptrCl->getAnalog(a); amplitude += a; nMMs++;
	}
	if (nMMs>=8) mMAvsP->Fill(amplitude/nMMs,fabs(1/H(5)));
      }
#endif
    }
    h_nBeams[0]->Fill((double)nBeams); h_nBoTs[0]->Fill((double)nBoTs);
    if      (evtTrig&stdTrigs) {
      h_nBeams[1]->Fill((double)nBeams); h_nBoTs[1]->Fill((double)nBoTs);
    }
    else if ((evtTrig&pureTrigs)==evtTrig) {
      h_nBeams[2]->Fill((double)nBeams); h_nBoTs[2]->Fill((double)nBoTs);
    }
    else if (evtTrig&specialTrigs) {
      h_nBeams[3]->Fill((double)nBeams); h_nBoTs[3]->Fill((double)nBoTs);
    }
#ifdef DEBUG_TRACKING
    nsp1t += nsp1; nsp2t += nsp2; nmust += nmus;
    if (TOpt::Print[8])
      printf("Evt %5d   %2d %2d %2d  %5d %5d %5d\n",
	     ievt,nsp1,nsp2,nmus,nsp1t,nsp2t,nmust);
#endif
  }

  if (TOpt::Hist[1]&0x10) {

    // ====================================================================
    // ********** RESIDUALS, EFFICIENCY, ``ON-TRACK'' DIGIT DATA **********
    // ...and else...
    // ====================================================================

    list<TTrack>::iterator it; int igr, ipl, idet;

    // For all of the above: one needs precise timing
    //  =>      ********** REQUIRE PRECISELY TIMED EVENTS **********

    const unsigned int allTrigs = 0xffff; // Although TCS can handle 0x7fffff, it was decided to limit max. #triggers to 16. And take advantage of the fact to make use of the other bits of the trigger TDC module.
    unsigned int evTrig = trig_mask&allTrigs; // Cut away trailing bits
    bool evtTCorr = false;
    // ***** CHECK WHETHER BADLY TIMED TRIGGER...
    // (For the check, we refer to the trigger described by ReMode[36], for
    // which a timing of the event by the beam track is required. This timing
    // may not be granted: if there are multiple beam tracks and the
    // processing of the event did not reach vertexing and hence the ambiguity
    // could not be raised.)
    unsigned int inclMask =// Trigger pattern for which bits we require only...
      TOpt::ReMode[36]&(~TOpt::ReMode[40]);// ...incl. match, cf. UpdateDrifts
    if ((evtTrig&TOpt::ReMode[36])==evtTrig ||
	(evtTrig&inclMask) ||
	// (Note: A null trigger pattern can happen in MC and will pass the
	// condition supra. It's probably not of any relevance here to address
	// the case at this point. Let's not do anything special to cope w/ it.)
	// ...Or alternative master trigger is defined (meaning also is needed).
	// Note also that in that case, "eventTime" is always, by definition,
	// defined (and hence, w/ a very high probability, finite, fulfilling
	// the 1st condition infra).
	ptrEvt()->getAlterMasterTrigger()>=0 ||
	// ...Or 2nd iteration of tracking going on
	ptrEvt()->reTrackingON()) {
      if (timingRequired &&
	  eventTime==0 &&                       // ... REQUIRE "eventTime" *****
	  eventTRef==0) return;                 //       ...or "eventTRef" *****
      if (eventTime) evtTCorr = true;
    }
    static bool imloBook = true;  // ***** PATTERNS of DET# ASSOCIATED TO IMLO
    static unsigned int imloMasks[4] = {0x1,0x102,0x4,0x608};
    static int imloHWords[4][2][2];       // 4 types of trigger x 2 hodos per...
    static unsigned imloHPats[4][2][2];   // trigger x 2 words for each hodo
    if (imloBook) {
      imloBook = false;
      int imlo, ih, iw; for (imlo = 0; imlo<4; imlo++) for (ih = 0; ih<2; ih++)
	for (iw = 0; iw<2; iw++) {
	  imloHWords[imlo][ih][iw] = -1; imloHPats[imlo][ih][iw] = 0;
	}
      if ((int)setup.vIplFirst().size()>=4) {
	for (int ipl = setup.vIplFirst()[2]; ipl<setup.vIplLast()[3]; ipl++) {
	  const TDetect &d = setup.iPlane2Detect(ipl);
	  if (d.Name[0]!='H') continue;
	  int w = ipl/32;
	  if (d.Name[1]=='O') {
	    imlo = 3; ih = d.Name[3]=='3' ? 0 : 1;  // HO03 and HO04
	  }
	  else {
	    ih = d.Name[3]=='4' ? 0 : 1;            // H?04 and HO05
	    if      (d.Name[1]=='I') imlo = 0;
	    else if (d.Name[1]=='M') imlo = 1;
	    else if (d.Name[1]=='L') imlo = 2;
	    else {
	      CsErrLog::msg(elError,__FILE__,__LINE__,
		"Detector \"%s\" is an hodo but none of IMLO",d.Name.c_str());
	      continue;
	    }
	  }
	  int jw; for (jw = 0, iw = -1; jw<2; jw++) {
	    if (w==imloHWords[imlo][ih][jw] ||
		imloHWords[imlo][ih][jw]==-1) { iw = jw; break; }
	  }
	  if (iw==-1)
	    CsErrLog::msg(elFatal,__FILE__,__LINE__,
   "Error assigning IMLO pattern to \"%s\" ipl=%d: words %d, %d already used",
	    d.Name.c_str(),ipl,imloHWords[imlo][ih][0],imloHWords[imlo][ih][1]);
	  imloHWords[imlo][ih][iw] = w; imloHPats[imlo][ih][iw] |= 1<<ipl%32;
	}
      }
    }
    int iMLO, imlo; for (imlo = 0, iMLO = -1; imlo<4; imlo++) {
      if (evtTrig&imloMasks[imlo]) { iMLO = imlo; break; }
    }
    THlx *loneBeam; for (it = listTrack.begin(), loneBeam = 0;
			 it!=listTrack.end(); it++) {
      TTrack &t = *it; int igr_t = t.Type;
      if (t.Type==0x10) {
	if (loneBeam) { loneBeam = 0; break; }
	else            loneBeam = &(t.Hlast);
      }
    }

#ifdef RDMon_GEMWrBack
    int GEMs_wrback = false;
#endif

    for (it = listTrack.begin(); it!=listTrack.end(); it++) {

      // **********************************************************************
      // ######################### LOOP OVER TRACKS #########################
      // **********************************************************************

      TTrack &t = *it; int igr_t = t.Type;
      double evtT = evtTCorr? eventTime : 0;
#if ! defined RDMon_DIGIT_DATA || RDMon_DIGIT_DATA < 2
      if (fabs(t.MeanTime-evtT-eventTRef)/t.SigmaTime>5./3)         // Off time?
	// Note: When timing undef: "SigmaTime"<0
	continue;
      // Else some particular studies do not care about timing => Allow all tracks , in order to maximize statistics
#endif
      if ((igr_t&TOpt::Hist[18])!=TOpt::Hist[18]) continue;
      if (TOpt::ReMode[3]==0 &&           // Track bridging is ON...
	  t.NGroups()<2 && 	          // ... => Require bridging...
	  TOpt::Hist[18]!=16) continue;   // ...unless zones==0x1
      if (!t.Hfirst.with_mom()) {
	CsErrLog::msg(elError,__FILE__,__LINE__,
		      "Run %d Evt %d Track #d spans %d groups but momentumless",
		      run,event,t.NGroups()); continue;
      }
      if (t.Type==0xc) continue;

      THlx Hstrt, Hextr;

      //     ********** ACCIDENTALLY COINCIDENT HALO TRACK? **********
      bool haloTrack = false;
      // At this point may remain accidentally coincident halo tracks,
      // travelling in detectors w/ poor time resolution (=> non negligible
      // probability of passing time cut supra)
      //  =>  I) Indentify those halo tracks.
      //     => Require a lone scifiSi track, to assign to the beam particle.
      //     II) They may still be usefull to probe the detectors downstream
      //        of muWall#2, on the Salève side: many of the tracks ending up
      //        there are halos triggering the outer trigger.
      //     => If some of the probed detectors lie downstream of muWall#2 (
      //       i.e. belong to zone 0x8): retain the halo tracks, while
      //       cutting on their sigmaTime.
      //        Else, throw away.
      //     N.B.: Still a possibility to throw them away in any case by
      //     defining "RDMon_DISCARD_HALOS"
      if (fabs(t.Hfirst.Mom())>TOpt::dCut[4]*.5) {// High P (=80 for 160 GeV)
	THlx &hFirst = t.Hfirst;
	if (fabs(hFirst(3))<.025 && fabs(hFirst(4))<.025) { // ...Paraxial...
	  //                         ***** CANDIDATE ACCIDENTAL HALO TRACK *****
	  int isMup = iMLO!=-1 && (t.Type&0x8);
	  if (isMup) {             // ***** IF EVENT IS IMLO: TRACK =? mu' *****
	    list<int>::const_iterator ipl, ihp; int hodos;
	    for (ipl = t.lPlnRef.begin(), ihp = t.lHitPat.begin(), hodos = 0;
		 ipl!=t.lPlnRef.end(); ipl++, ihp++) {
	      if (*ihp<0) continue;
	      int jpl = *ipl, w = jpl/32, bit = 1<<jpl%32;
	      for (int ih = 0; ih<2; ih++) for (int iw = 0; iw<2; iw++) {
		if (w!=imloHWords[iMLO][ih][iw]) continue;
		if (bit&imloHPats[iMLO][ih][iw]) hodos |= 1<<ih;
	      }
	    }
	    if (hodos!=3) isMup = false;  
	  }
	  if (!isMup) {                                // ***** IF !mu'... *****
	    if (loneBeam) {                     // ...REQUIRE LONE BEAM... *****
	      Hextr = t.Hfirst; THlx Hbeam = *loneBeam;
	      double X0 = setup.TargetCenter[0];
	      if (Hextr.FindCDA(Hbeam,X0,X0-200,X0+200)) {  // ...FIND CDA *****
		double dy = Hextr(1)-Hbeam(1), dz = Hextr(2)-Hbeam(2);
		double r = sqrt(dy*dy+dz*dz);
		double dr = sqrt((Hextr(1,1)+Hbeam(1,1))*dy*dy+
				 (Hextr(2,2)+Hbeam(2,2))*dz*dz+
				 2*(Hextr(1,2)-Hbeam(1,2))*dy*dz)/r;
		if (r/dr>2) haloTrack = true;          // ...CUT @ 2*sigma *****
	      }
	    }
	  }
	}
      }  // End indentifying halo track
#ifdef RDMon_DISCARD_HALOS         // ***** DISCARD HALO TRACK UPON OPTION *****
      if (haloTrack) continue;
#endif

#ifdef DEBUG
      int debug = 0;
#endif
#ifdef RDMon_AligTree
      nhAlig = 0; int iAlig_na[2] = {0,0};
#endif
#ifdef RDMon_SM2Tree
      static int sm2_ipl; sm2_pguess = 0; sm2_pfound = 0;
#endif
#ifdef RDMon_GEMWrBack
      int nhGEMs_2 = 0, nhAll_2 = 0;
#endif

      // VS,S,LAT: iterators on hits and planes of TTrack
      list<int>::const_iterator ip_t = t.lPRef().begin(); 
      list<int>::const_iterator ih_t = t.lHPat().begin(); 
      int vsslatz = 0; // Init VS,S,LAT*zones description of TTrack

      for (igr = 0, idet = 0; igr<igr_mx; igr++) {
	if ((1<<igr&(TOpt::Hist[17]|    // Plane to monitor if inactive
		     TOpt::Hist[16]))   // Plane to monitor whatsoever
	    ==0) continue;
	int resId = (1<<igr)&igr_t ? 0 : 1; // Type of residuals to be filled
	if (resId==1 &&
	    (igr==3 && (igr_t&0x4)==0 || // Don't extrap. SM1 track to zone 0x4
	     igr==4 && (igr_t&0x1)==0))  // Don't extrap. SM2 track to zone 0x8
	  continue;

	// ***** VECTOR of HELICES to be USED for this->ZONE *****
	set<THlx0,less<THlx0> > hlx0s;
	double Xfirst = setup.iPlane2Detect(setup.vIplFirst()[igr]).X(0)-.1;
	double Xlast  = setup.iPlane2Detect(setup.vIplLast()[igr]) .X(0)+.1;
	if (igr==2 && TOpt::MuonWall[10]>Xlast) Xlast = TOpt::MuonWall[10];
	if (Xfirst<=t.Hfirst(0) && t.Hfirst(0)<=Xlast)
	  hlx0s.insert(THlx0(t.Hfirst));                 // 1ST HELIX?
	if (Xfirst<=t.Hlast(0) && t.Hlast(0)<=Xlast)
	  hlx0s.insert(THlx0(t.Hlast));                  // LAST HELIX?
	list<THlx>::const_iterator iHsm;                    // SMOOTHED HELICES?
	for (iHsm = t.lHsm.begin(); iHsm!=t.lHsm.end(); iHsm++) {
	  if (Xfirst<=(*iHsm)(0) && (*iHsm)(0)<=Xlast)
	    hlx0s.insert(THlx0(*iHsm));
	}
	if (hlx0s.empty() && resId==1) {
	  // Rescue Hfirst/Hlast when extrapolating outside track's range
	  if      (t.Hfirst(0)>Xlast) hlx0s.insert(THlx0(t.Hfirst));
	  else if (Xfirst>t.Hlast(0)) hlx0s.insert(THlx0(t.Hlast));
	}
	if (hlx0s.empty()) {
	  // This means track includes zone where dets are to be monitored
	  // And hence that it could have been smoothed.
	  // => Require that a smoothed point be booked!
	  CsErrLog::
	    msg(elError,__FILE__,__LINE__,
		"TTrack #%d 0x%x has no helix (1st,last,smooth) for zone 0x%x",
		t.Id,igr_t,1<<igr);
	  continue;
	}
	//	hlx0s.sort();
	Hstrt = (*hlx0s.begin()).hlx;
	
	set<THlx0,less<THlx0> >::const_iterator currentHlx0 = hlx0s.begin(),
	  nextHlx0 = ++currentHlx0;
	double currentHlX0 = Hstrt(0), nextHlX0 =
	  nextHlx0!=hlx0s.end() ? (*nextHlx0).hlx(0) : 1.e6;
	//#define RDMon_BEAM_WITH_MOM
#ifdef RDMon_BEAM_WITH_MOM
	if (!Hstrt.with_mom() &&
	    (igr==4 && igr_t==0x10 ||   // Extrapolating thru target solenoid
	     igr==0 && igr_t==0x1 /* e.g. from single zone reco */))
	  Hstrt(5) = dCut[4];  // WARNING: Goes wrong if next helix is used!
#endif

	// ############### LOOP OVER ``BOOKED'' GROUPS ###############

	if (igr<4 && TOpt::Hist[16]&0x20) {

	  // ##### EXTRAPOLATIONS TO ORIGINS #####

	  int iintrcpt = -1;
	  if   (igr==0) {                  // Zone 0..
	    Hextr(0) = setup.TargetCenter[0];    // Origin is target
	    iintrcpt = (igr_t&0x3)==0x3 ? 0 : 1; // Histo[0] = 1st spectro
	  }
	  else if (igr==1 && (igr_t&0x3)==0x3 ||
		   igr==2 && (igr_t&0x6)==0x6 ||
		   igr==3 && igr_t&0x8) {
	    Hextr(0) =                     // Else Origin is beginning of zone
	      setup.iPlane2Detect(setup.vIplFirst()[igr]).X(0);
	    iintrcpt = igr+1;                 // Histo[2] = zone 1
	  }
	  if (iintrcpt>=0 && Hstrt.Extrapolate(Hextr)) {
	    if      (igr<3) intrcpts[iintrcpt]->Fill(Hextr(1),Hextr(2));
	    else { //if (igr==3)
	      if (t.radLenFr>80) intrcpts[4]->Fill(Hextr(1),Hextr(2));
	      else {
		intrcpts[5]->Fill(Hextr(1),Hextr(2));
		if (fabs(Hextr(1)-30)>30 || fabs(Hextr(2))>30) {
		  printf("Id%d 0x%x %f %f XX0 %f\n",t.Id,t.Type,
			 Hextr(1),Hextr(2),t.radLenFr);
		}
	      }
	    }
	  }
	}


#ifdef RDMon_DBL_LAYER
	float dblxprv = 0; double dblr[2] = {0,0}, dblrprv[2] = {0,0};
	static float Uguessprv, Vguessprv, Urndomprv;
	static int firedprv, accidentprv;
#endif
#ifdef RDMon_DIGIT_GM
	float gemxprv = 0; double gemaprv = 0;
#endif
	//#define RDMon_EXTRAP_BEAM
#ifdef RDMon_EXTRAP_BEAM
	// Aimed at debugging momentum assignment to tracks extrapolated
	// through target
	bool first_4 = true;
#endif
	//#define RDMon_REFIT_HALO
#ifdef RDMon_REFIT_HALO
	int n_halo = 0, id_halo = idet;
#endif
#ifdef RDmon_SCIFI_TIME
	double scifitimes[7]; int n_scifi = 0;
#endif
	for (ipl = setup.vIplFirst()[igr];
	     ipl<=setup.vIplLast()[igr]; ipl++) {

	  // ************************************************************
	  // #################### LOOP OVER PLANES ####################
	  // ************************************************************

	  const TPlane &p = setup.vPlane(ipl);
	  const TDetect &d = setup.vDetect(p.IDetRef);

	  while (ip_t!=t.lPRef().end() && *ip_t<ipl) {
	    ip_t++; ih_t++;    // Update ``track'' iterators
	  }
	  if (igr<3) { // Excluding zones 3: beyond wall, and 4: before target
	    // VS,S,LAT: Update description...
	    if (*ip_t==ipl && *ih_t>=0) {  // ... if plane belongs to track
	      if      (d.IType==22 || d.IType==21)
		/*       */         vsslatz |= 0x1<<igr*3; // SciFi/Si => VSAT
	      else if (d.IType==26 ||                      // GEMs or ..
		       d.IType==27) vsslatz |= 0x2<<igr*3; // ...uMs => SAT
	      else                  vsslatz |= 0x4<<igr*3; // Otherwise LAT
#ifdef RDMon_GEMWrBack
	      if (igr==2) {
		if (d.IType!=26) nhGEMs_2++;
		nhAll_2++;
	      }
#endif
	    }
	  }

	  // ********** IS CURRENT PLANE TO BE MONITORED? **********
#ifdef RDMon_DIGIT_GM                            // Special case of GEM/Si amps
	  if (d.IType!=21 && d.IType!=26 && d.IType!=28)
	    continue;
#else
	  if (p.IFlag && (1<<igr&TOpt::Hist[16])==0 &&  // ... or Det's !off
	      (igr!=4 || igr_t!=16))            // Special case of beam tracks
	    continue;
#endif
	  
	  //  ***** REQUIRE TRACK TO FIRE SOME of CURRENT TStation PLANES *****
	  const TStation *&st = p.Station;
	  int spli = st->IPlanes[0], splf = spli+(int)st->IPlanes.size()-1;
	  int nFiredPlanes = 0;     // Count # of firing planes in TStation
	  if (*ip_t==ipl) {
	    list<int>::const_iterator jp_t, jh_t, kp_t = t.lPRef().begin();
	    jp_t = ip_t; jh_t = ih_t; // Start from current plane, descending
	    do {
	      if (*jp_t!=ipl &&// Skip current plane (useful in the GEM/Si case)
		  *jh_t>=0) nFiredPlanes++;
	      if (jp_t!=kp_t) { jp_t--; jh_t--; }
	    } while (jp_t!=kp_t && *jp_t>=spli);
	    jp_t = ip_t; jh_t = ih_t; jp_t++; jh_t++; // And now ascending
	    while (jp_t!=t.lPRef().end() && *jp_t<=splf) {
	      if (*jh_t>=0) nFiredPlanes++;
	      jp_t++; jh_t++;
	    }
	  }
	  if (nFiredPlanes==0 && d.IType<41) continue; // EXCLUDING HODOS (Note: would be good to exclude also GM11 (others?)
#if defined RDMon_DEBUG_MAs && RDMon_DEBUG_MAs > 1 // Dump event w/ track in MAs
	  if (d.IType==39 || d.IType==40) {
	    static int iStream = -1; if (iStream<0)
	      iStream = CsEvent::Instance()->openRawEventOutputStream("ManyMAs.raw");
	    CsEvent::Instance()->outputRawEventToStream(iStream);
#  if RDMon_DEBUG_MAs > 2
	    static int nMAEvts = 0;
	    printf("RDMon_DEBUG_MAS Evt %d Trk %d  %d\n",
		   ptrEvt()->getEventNumberInRun(),t.Id,++nMAEvts);
#  endif
	  }
#endif


	  double currentX = d.X(0);
	  while (fabs(currentX-nextHlX0)<fabs(currentX-currentHlX0)) {
	    // ***** SWITCH ON NEXT AVAILABLE HELIX IF CLOSER *****
	    Hstrt = (*nextHlx0).hlx;
	    currentHlx0 = nextHlx0; nextHlx0 = ++currentHlx0;
	    currentHlX0 = Hstrt(0);
	    nextHlX0 = nextHlx0!=hlx0s.end() ? (*nextHlx0).hlx(0) : 1.e6;
	  }

	  Hextr(0) = currentX; // ********** EXTRAP. TO CURRENT PLANE **********
	  if (!Hstrt.Extrapolate(Hextr)) { idet++; continue; }// Extrapolate
	  Hstrt = Hextr;
#ifdef RDMon_EXTRAP_BEAM
	  if (igr==4 && first_4) { t.Hfirst = Hextr; first_4 = false; }
#endif

	  // ********** CHECK in ACTIVE **********
	  if (!d.InActive(Hextr(1),Hextr(2))) { idet++; continue; }

	  // ********** CHECK in TIME for CURRENT PLANE TYPE **********
	  int isDrift = 0;
	  if (d.PtrDet()->hasDrift()) {  // Guarantee that not in MWPC mode
	    // Stricter than PtrDet()->hasDrift? Has to be set in agrement
	    // with ``drift'' definition in the booking block supra
	    if (strncmp(d.Name.c_str(),"DC",2)==0) isDrift = 1;
	    if (strncmp(d.Name.c_str(),"ST",2)==0) isDrift = 2;
	    if (strncmp(d.Name.c_str(),"DW",2)==0) isDrift = 3;
	    if (strncmp(d.Name.c_str(),"MB",2)==0) isDrift = 4;
	  }
	  if (timingRequired && 
	      fabs(t.SigmaTime)>trackTMx[isDrift]) { idet++; continue; }

	  // ********** INTERSECTION w/ CURRENT PLANE ...**********
	  double Uguess = d.Ca*Hextr(1)+d.Sa*Hextr(2);
	  double deltaU = sqrt(d.Ca*d.Ca*Hextr(1,1)+
			       d.Sa*d.Sa*Hextr(2,2)+
			     2*d.Ca*d.Sa*Hextr(1,2));
	  double Vguess = d.Ca*Hextr(2)-d.Sa*Hextr(1);
	  // ... and RANDOM
	  double Urndom =  d.X(1)*d.Ca+d.X(2)*d.Sa +
	    d.Siz(1)*(CsRandom::flat()-.5);

	  //   *************** SEARCH for MATCHING CLUSTERS ***************
	  double du_best, du_rndom, du_good, du_effic; int n_good = 0;
	  static THit *h_best, *h_rndom; int fired, accident;
	  // "du_good": Cut loose enough to have retained hit span a large
	  //           enough domain about residuals peak so that one can
	  //           make out the flat background of random coinc.
#ifdef RDMon_DEBUG_MAs
	  // While debuggging MAs (i.e. trying to determine the extent and
	  // position of their dead zone), we are not interested in doing a
	  // real efficiency measurement w/ noise subtraction => let's
	  // restrict as much as possible the search route.
	  du_good = d.Resol*6;
#else
	  du_good = d.Resol*10;
#endif
	  // "du_effic": Tighter cut, used for efficiency determination. A
	  //            tight enough cut is needed to get the measured efficieny
	  //            to account for possible misalignment, even though the
	  //            the contribution of random coinc, that can become
	  //            significant w/ the somewhat loose "du_good" cut, can be
	  //            corrected for using the measured background rate.
	  //             What value should one take, is not staightforward to
	  //            tell. In principle, a value equal to the width of the
	  //            route, expressed in multiple of the detector's
	  //            resolution, opened in the track finding algorithm.
	  //            Problem: this is not unambiguously defined: it depends
	  //            upon the detector's zone, and, above all, varies in the
	  //            course of the algorithm, with increasing iteration #.
	  du_effic = d.Resol*3;
#ifdef VERY_LOW_RESOLUTION
	  du_good *= 4;
#endif
	  du_best=du_rndom = du_good;
      
#ifdef RDMon_SM2Tree
	  if (igr==0) {
	    sm2_ipl = ipl-setup.vIplFirst()[0];
	    if (sm2_ipl<12) {
	      sm2_uguess[sm2_ipl] = Uguess; sm2_vguess[sm2_ipl] = Vguess;
	      sm2_pguess |= (1<<sm2_ipl);
	    }
	    else CsErrLog::mes(elFatal,"RDMon_SM2Tree: Too many dets!");
	  }
#endif

	  vector<int>::const_iterator i;
	  for (i = p.vHitRef().begin(), fired=accident = 0;
	       i!=p.vHitRef().end(); i++) {

	    // ********** LOOP OVER HIT REFERENCES ON THE PLANE **********

	    THit &h = vecHit[*i]; // ref to Thit vector element

	    if (isDrift && p.IFlag==0 &&  // If Drift ON: correction yet done
		TOpt::ReMode[29]&0x2 && !this->isMC) {

	      // ***** CORRECT PROPAGATION in DRIFTS *****

	      double tT = 0;
	      if (haloTrack)     tT = t.MeanTime;
	      else if (evtTCorr) tT = eventTime; 
	      double y = Hextr(1), z = Hextr(2);
	      CsDetector *det = d.PtrDet();
	      bool error; double U = det->getCorrU(h.PtrClus(),y*10,z*10,tT,error);
	      if (!error) h.u = U/10;
	      else continue;
	    }        // End correct for propagation in drifts

#define PLOT_ALL_RESIDUALS
#ifdef PLOT_ALL_RESIDUALS
	    // ********** PLOT ALL HITS **********

	    if (isDrift && p.IFlag) {
	      if (h.ptrCl->getLRProb()>.5) {
		// If OFF AND out of LR: proba IS 50% and infra to be  skipped
		// => histo is filled later on
		((residus[resId])[idet])->Fill(Uguess,h.U-Uguess);
		((reperps[resId])[idet])->Fill(Vguess,h.U-Uguess);
	      }
	    }
	    else {
	      ((residus[resId])[idet])->Fill(Uguess,h.U-Uguess);
	      ((reperps[resId])[idet])->Fill(Vguess,h.U-Uguess);
	    }
#endif

	    // ********** BEST HIT **********

	    if (fabs(h.U-Uguess)<sqrt(du_good*du_good+9*deltaU*deltaU)) {
	      if (!n_good || fabs(h.U-Uguess)<du_best) {
		du_best = fabs(h.U-Uguess); h_best = &h;
		if (!fired) fired = 1;
		if (du_best<sqrt(du_effic*du_effic+9*deltaU*deltaU)) fired = 2;
	      }
	      n_good++;
	    }
	    if (fabs(h.U-Urndom)<sqrt(du_rndom*du_rndom+9*deltaU*deltaU)) {
	      du_rndom = fabs(h.U-Urndom); h_rndom = &h;
	      if (!accident) accident = 1;
	      if (du_rndom<sqrt(du_effic*du_effic+9*deltaU*deltaU))
		accident = 2;
	    }

	  }// *************** End of loop over hits on the plane ***************
#ifdef DEBUG
	  if (debug && fired)
	    cout<<d.Name<<" "<<Uguess<<" "<<h_best->U<<" "<<h_best->U-Uguess<<endl;
#endif

#ifdef RDMon_AligTree
	  if ((d.IType==26 || d.IType==27) &&  // Alig tree assumed GM=MM
	      n_good==1) {        	         // GM or MM AMBIGUITY?
	    int ialig = Alig_map[ipl-Alig_ipl0];
	    if (ialig>=0) {
	      iAlig_na[ialig/32] |= 1<<ialig%32; // Alig non-ambiguous pattern
	      nhAlig++; 
	    }
	  }
#endif

	  if (fired
#ifdef PLOT_ALL_RESIDUALS
	      // In case of drift, if it's off: probably means
	      // that LR has not been evaluated and no hit has passed
	      // the condition set on LR probability supra. Since one is
	      // probably interested in characterising detector rather
	      // than quality of LR, do as if the latter is ideal and
	      // plot the best (i.e. closest) hit
	      && isDrift && p.IFlag==0
#endif
	      ) {               // ********** PLOT BEST HIT **********
	    ((residus[resId])[idet])->Fill(Uguess,h_best->U-Uguess);
	    ((reperps[resId])[idet])->Fill(Vguess,h_best->U-Uguess);
	  }

#ifdef RDMon_SM2Tree
	  if (fired) {
	    sm2_found[sm2_ipl] = h_best->U;
	    sm2_pfound |= (1<<sm2_ipl);
	  }
#endif
#ifdef RDMon_REFIT_HALO
	  if (igr==4 && fired>1) { // Pick-up hits in 1st 2 detectors and refit
	    if (d.Name.find("FI")==0) {
	      t.AddHit(*h_best); n_halo ++;
	    }
	  }
#endif

	  // *************** EFFICIENCY DETERMINATION ***************
#ifdef DEBUG
	  int dump = 0; if (dump) t.DumpHits(dump-1);
#endif
	  if (fired>1) {
	    // Fill both "guess" and "effic" w/ "h_best", whereas "guess"
	    // will be filled w/ "Uguess" when det is inefficient. This is
	    // prefered over filling everybody w/ "guess", for it allows to
	    // better single out bad wires.
	    guess[idet]->Fill(h_best->U); effic[idet]->Fill(h_best->U);
#if defined RDMon_EFF_MAPs && !defined RDMon_DBL_LAYER
	    guess2[idet]->Fill(h_best->U,Vguess);
	    effic2[idet]->Fill(h_best->U,Vguess);
#endif
	  }
	  else {
	    guess[idet]->Fill(Uguess);
#if defined RDMon_EFF_MAPs && !defined RDMon_DBL_LAYER
	    guess2[idet]->Fill(Uguess,Vguess);
#endif
	  }
	  if (accident>1) {	          // ********** BACKGROUND **********
	    rndom[idet]->Fill(h_rndom->U); bckgd[idet]->Fill(h_rndom->U);
	  }
	  else rndom[idet]->Fill(Urndom);

#ifdef RDMon_DBL_LAYER
	  if (isDrift) {// ********** DBL LAYER RESOLUTION/EFFICIENCY **********
	    int hist = dblhists[d.IDet];
	    bool ok = true;
	      //#define RDMon_MB_REJECT_PARITY 3  // 1 <-> event, 2<-> odd, 3 <-> debug
#  ifdef RDMon_MB_REJECT_PARITY
	    static int channelNum;// Enlarged scope: only for debugging purposes
#  endif
	    if (fired) {
#  if defined RDMon_TOP_SPLIT_STRAWS || defined RDMon_MB_REJECT_PARITY
	      // Consider split straws separately. Whether it's top (most
	      // interesting, for we suspect T0 to be wrongly defined there) or
	      // everything but top that we want to consider is determined from
	      // RDMon_TOP_SPLIT_STRAWS
	      bool getDigit = false;
#    if   defined RDMon_TOP_SPLIT_STRAWS
	     if (strncmp(d.Name.c_str(),"ST",2)==0) getDigit = true;
#    elif defined RDMon_MB_REJECT_PARITY
	     if (strncmp(d.Name.c_str(),"MB",2)==0) getDigit = true;
#    endif
	     if (getDigit) {
		list<CsDigit*> lDig = h_best->ptrCl->getDigitsList();
		list<CsDigit*> ::iterator id;
		for (id = lDig.begin(); id!=lDig.end(); id++) {
		  // Loop over CsDigits in the CsCluster. Should be only one.
#    if   RDMon_TOP_SPLIT_STRAWS == 1     // Select top split (or is it bottom?)
		  if ((*id)->getData()[1]!=-1) ok = false;
#    elif RDMon_TOP_SPLIT_STRAWS == 2     // All but top split (or but buttom?)
		  if ((*id)->getData()[1]==-1) ok = false;
#    elif RDMon_TOP_SPLIT_STRAWS == 3     // Bottom split (or top?)
		  if ((*id)->getData()[1]!=1) ok = false;
#    elif defined RDMon_MB_REJECT_PARITY
		  channelNum = (*id)->getAddress();
		  if (channelNum%2+1==RDMon_MB_REJECT_PARITY) ok = false;
#    endif
		}
	      }
#  endif
	      // For all what follows: do not care about LR probabilities
	      // for it's rather a property of tracking while we are
	      // trying to characterise detectors. Instead, keep track
	      // of both genuine and mirror.
	      dblr[0] = h_best->U-Uguess;
	      int mirror; if ((mirror = h_best->Mirror)!=-1) {
		dblr[1] = vecHit[mirror].U-Uguess;
	      }
	      else dblr[1] = h_best->U-Uguess;
	    }
	    double dx_mx = 1.5;  // Max distance between 2 layers of DL
	    if      (d.Name[0]=='M') dx_mx = 3;
	    else if (d.Name[1]=='W') dx_mx = 2.5;
	    if (fabs(d.X(0)-dblxprv)<dx_mx) {      // 2nd layer: Fill
	      int igm, jgm; double diff;
	      for (igm = 0, diff = 100; igm<2; igm++)
		for (jgm = 0; jgm<2; jgm++) {
		  double d = dblr[igm]-dblrprv[jgm];
		  if (firedprv && fired && ok)
		    dblreso2[hist/2]->Fill(d,h_best->U);
		  if (fabs(d)<fabs(diff)) diff = d;
		}
	      if (firedprv && fired && ok)
		dblresol[hist/2]->Fill(diff,h_best->U);
	      if (firedprv>1 || fired>1) {
		double bestU = 0;
		if (fired>1 && firedprv>1)
		  bestU = (dblr[0]+dblrprv[0]+Uguess+Uguessprv)/2;
		else
		  bestU = fired ? dblr[0]+Uguess : dblrprv[0]+Uguessprv;
		dbleffic[hist/2]->Fill(bestU);
		dblguess[hist/2]->Fill(bestU);
#ifdef RDMon_EFF_MAPs
		guess2[idet]->Fill(bestU,(Vguess+Vguessprv)/2);
		effic2[idet]->Fill(bestU,(Vguess+Vguessprv)/2);
#endif
	      }
	      else {
		dblguess[hist/2]->Fill((Uguess+Uguessprv)/2);
#ifdef RDMon_EFF_MAPs
		guess2[idet]->Fill((Uguess+Uguessprv)/2,(Vguess+Vguessprv)/2);
#endif
	      }
	      if (accident>1 || accidentprv>1) {
		bool prv = accidentprv>1 ? true : false;
		if (accident>1 && accidentprv>1)   // If both: only one entry
		  if (CsRandom::flat()>.5) prv = false;
		if (prv) {
		  dblbckgd[hist/2]->Fill(Urndomprv);
		  dblrndom[hist/2]->Fill(Urndomprv);
		}
		else {
		  dblbckgd[hist/2]->Fill(Urndom);
		  dblrndom[hist/2]->Fill(Urndom);
		}
	      }
	      else {
		if (CsRandom::flat()>.5) dblrndom[hist/2]->Fill(Urndom);
		else                     dblrndom[hist/2]->Fill(Urndomprv);
	      }
	    }
	    else {                              // 1st layer: Save résidus
	      dblrprv[0] = dblr[0]; dblrprv[1] = dblr[1];
	      Uguessprv = Uguess; Vguessprv = Vguess; Urndomprv = Urndom;
	      firedprv = fired; accidentprv = accident;
	    }
	    dblxprv = d.X(0);           // Update "dblxprv"
	  }
#endif

	  // ********** DIGIT (or CLUSTER) DATA **********

#ifdef RDMon_DIGIT_DATA
	  //                               ***** FILL HISTOS OF DIGIT DATA *****
#  ifdef RDMon_DIGIT_MM
	  if (strncmp(d.Name.c_str(),"MM",2)==0) {
	    if (fired) {
	      int hist = mmhists[d.IDet];
	      double ct; h_best->ptrCl->getTime(ct);
	      mmctime[hist]->Fill(ct);
	      double cToT; h_best->ptrCl->getAnalog(cToT);
	      mmcToT[hist]->Fill(cToT);
	      list<CsDigit*> lDig = h_best->ptrCl->getDigitsList();
	      list<CsDigit*> ::iterator id;
	      for (id = lDig.begin(); id!=lDig.end(); id++) {
		// loop over CsDigits in the CsCluster
		static const float leadtWght = 0.6;      
		double ht = (*id)->getData()[0];
		int wire = (*id)->getAddress();
		mmhtime[hist]->Fill((double)wire,ht);
		double hToT  = (*id)->getData()[1];
		double lt = ht-(1-leadtWght)*hToT;
		mmltime[hist]->Fill(lt);
		double tt = ht+leadtWght*hToT;
		mmttime[hist]->Fill(tt);
		if (384<=wire && wire<640)
		  mmhtvT[hist]->Fill(hToT,ht);
	      }
	    }
	  }
#  endif
#  ifdef RDMon_DIGIT_DC
	  if (isDrift) {
	    int hist = drifthists[d.IDet];
	    int icell = (int)((Uguess-d.Uorig-d.Pitch/2)/d.Pitch);
	    double R = Uguess-d.Uorig-d.Pitch-icell*d.Pitch;
	    if (fired>1) {
	      list<CsDigit*> lDig = h_best->ptrCl->getDigitsList();
	      list<CsDigit*> ::iterator id;
	      for (id = lDig.begin(); id!=lDig.end(); id++) {
		// Loop over CsDigits in the CsCluster
#  ifdef RDMon_TOP_SPLIT_STRAWS
		if (strncmp(d.Name.c_str(),"ST",2)==0) {
		  // Consider split straws separately.
#    if   RDMon_TOP_SPLIT_STRAWS == 1     // Select top split (or is it bottom?)
		  if ((*id)->getData()[1]!=-1) continue;
#    elif RDMon_TOP_SPLIT_STRAWS == 2     // All but top split (or but buttom?)
		  if ((*id)->getData()[1]==-1) continue;
#    elif RDMon_TOP_SPLIT_STRAWS == 3     // Bottom split (or top?)
		  if ((*id)->getData()[1]!=1) continue;
#    endif
		}
#  endif
		double T = (*id)->getDatum(); drifttime[hist]->Fill(T,R);
	      }
	      driftres[hist]->Fill(R,h_best->U-Uguess);
	      drifteff[hist]->Fill(R); 
	    }
	    driftges[hist]->Fill(R);
	    icell = (int)((Urndom-d.Uorig-d.Pitch/2)/d.Pitch);
	    R = Urndom-d.Uorig-d.Pitch-icell*d.Pitch;
	    if (accident>1) driftbck[hist]->Fill(R); 
	    driftrdm[hist]->Fill(R);
	  }
#  endif
#  ifdef RDMon_DIGIT_GM
	  if (strncmp(d.Name.c_str(),"GM",2)==0 ||
	      strncmp(d.Name.c_str(),"GP",2)==0 && d.Name[4]!='P' ||
	      strncmp(d.Name.c_str(),"SI",2)==0) {
	    int hist = gemhists[d.IDet];
	    if (fired>1) {
	      // HIT amplitude distributions (there are 3)
	      list<CsDigit*> lDig = h_best->ptrCl->getDigitsList();
	      list<CsDigit*> ::iterator id;
	      static int best_wire; double best_ampl3;
	      for (id = lDig.begin(), best_ampl3 = 0; id!=lDig.end(); id++) {
		// Loop over CsDigits in the CsCluster...
		int wire = (*id)->getAddress();
		double ampl3 = (*id)->getData()[2];
		if (ampl3>best_ampl3) {
		  // ...looking for wire with largest 3rd ampltude
		  best_ampl3 = ampl3; best_wire = wire;
		}
	      }
	      for (id = lDig.begin(); id!=lDig.end(); id++) {
		// Loop over CsDigits in the CsCluster...
		int wire = (*id)->getAddress();
		if (abs(wire-best_wire)<1) {
		  // ...and fill amplitudes for wire == largest +/- 0
		  double *ampls = (*id)->getData();
		  for (int iampl = 0; iampl<3; iampl++) {
#    ifdef RDMon_GM_2D_CALIB
		    (gemampls[iampl])[hist]->Fill((double)wire,ampls[iampl]);
#    else
		    (gemampls[iampl])[hist]->Fill(ampls[iampl]);
#    endif
		  }
		  gemratio[hist]->Fill(ampls[1]/ampls[2],ampls[0]/ampls[2]);
		}
	      }
	      // CLUSTER amplitude correlation
	      double amplc = h_best->ptrCl->getAllAnalogData()[2];
	      if (fabs(d.X(0)-gemxprv)<.1) {
		amplc = gemcorrs[hist/2][0] + gemcorrs[hist/2][1]*amplc+
		   gemcorrs[hist/2][2]*amplc*amplc;
		gemacorr[hist/2]->Fill(gemaprv,amplc);
		gemadiff[hist/2]->Fill(amplc-gemaprv);
	      }
	      else gemaprv = amplc;
	      double time; h_best->ptrCl->getTime(time);
	      gemctime[hist]->Fill(time);
	      //#  define RDMon_DEBUG_APV
#  ifdef RDMon_DEBUG_APV
	      vector<double> analogs = h_best->ptrCl->getAllAnalogData();
	      double si1 = analogs[1], si2 = analogs[2], si3 = analogs[3];
	      double tcs = CsEvent::Instance()->getTCSPhaseTime();
	      double tt  = CS::ChipF1::GetTT().GetTimeNorm();
#  endif
	      gemxprv = d.X(0);   // Update "gemxprv"
	    } // End of GM DIGIT_DATA if fired
	    else {
	      if (fabs(d.X(0)-gemxprv)<.1) {   // GEM did not fire...
		// ... but counterpart just encountered may have
		gemacorr[hist/2]->Fill(gemaprv,0.);
	      }
	      else // ... or counterpart still ahead could
		gemaprv = 0.;
	      gemxprv = d.X(0);   // Update "gemxprv"
	    }
	  }
#  endif
#endif
#ifdef RDmon_SCIFI_TIME
	  if (fired>1 && d.IType==22) {
	    double time; h_best->ptrCl->getTime(time);
	    scifitimes[n_scifi++] = time;
	  }
#endif
	  idet++;
	}        // ########## End loop on planes ##########

#ifdef RDmon_SCIFI_TIME
	if (n_scifi>=4) {
	  for (int i = 0; i<n_scifi; i++)
	    for (int j = i+1; j<n_scifi; j++) {
	      scifihists[igr]->Fill(scifitimes[j]-scifitimes[i]);
	    }
	}
#endif
#ifdef RDMon_REFIT_HALO
	if (n_halo>=5) {
	  int iter, ret;
	  for (iter = 0; iter<2; iter++) {
	    ret = iter==0 ? t.QuickKF(-1,1) : t.FullKF(-1);// backward
	    if(!ret) {
	      CsErrLog::mes(elWarning,
			    "Backward Kalman fit failed => No Bridging");
	      break;
	    }
	    ret = iter==0 ? t.QuickKF( 1,1) : t.FullKF( 1); // forward
	    if (!ret) {
	      CsErrLog::mes(elWarning,
			    "Forward Kalman fit failed => No Bridging");
	      break;
	    }
	  }
	  if (ret) {   // Extrapolate to SI and Fill residuals[0]
	    Hstrt = t.Hfirst;
	    for (ipl = setup.vIplFirst()[4], idet = id_halo;
		 ipl<=setup.vIplLast()[4]; ipl++, idet++) {
	      // ########## LOOP OVER PLANES ##########
	      const TPlane &p = setup.vPlane(ipl);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      if (d.Name.find("SI")) continue;
	      Hextr(0) = d.X(0);
	      if (!Hstrt.Extrapolate(Hextr)) continue;// Extrapolate
	      Hstrt = Hextr;
	      // ********** INTERSECTION w/ CURRENT PLANE ...**********
	      double Uguess = d.Ca*Hextr(1)+d.Sa*Hextr(2);
	      double Vguess = d.Ca*Hextr(2)-d.Sa*Hextr(1);
	      vector<int>::const_iterator i;
	      for (i = p.vHitRef().begin(); i!=p.vHitRef().end(); i++){
		// ********** LOOP OVER HIT REFERENCES ON THE PLANE **********
		THit &h = vecHit[*i]; // ref to Thit vector element
		// ***** PLOT ALL HITS *****
		((residus[0])[idet])->Fill(Uguess,h.U-Uguess);
		((reperps[0])[idet])->Fill(Vguess,h.U-Uguess);
	      }  // End loop on SI hits
	    }  // End loop on planes in search of SIs
	  }  // End case spectro track extended to FIs ante target fits
	}  // End spectro track extended to at least 5 FIs
#endif

      }     // ############### End loop on groups ###############

#ifdef RDMon_AligTree
    {
      nhAlig=nAlig[0]=nAlig[1]=nAlig[2]=nAlig[3] = 0; iAlig[0]=iAlig[1] = 0;
      list<int>::const_iterator ih = t.lHPat().begin();
      while (ih!=t.lHPat().end()) {
	if (*ih>=0) {
	  THit &h = vecHit[*ih];
	  int ialig = Alig_map[h.IPlane-Alig_ipl0];
	  if (ialig>=0) {
	    if (iAlig_na[ialig/32]&(1<<ialig%32)) {
	      iAlig[ialig/32] |= 1<<ialig%32;     // GEM hit pattern
	      uAlig[nhAlig] = h.U;                // Fill aray of Alig hits 
	      nhAlig++;                           // Count Alig hits
	      const TPlane &p = setup.vPlane(h.IPlane);
	      const TDetect &d = setup.vDetect(p.IDetRef);
	      if      (d.X(0)< 350) nAlig[0]++;
	      else if (d.X(0)< 750) nAlig[1]++;
	      else if (d.X(0)<1750) nAlig[2]++;
	      else                  nAlig[3]++;
	    }
	  }
	}
	ih++;
      }
      //#  define DEBUG_RDMon_AligTree
#  ifdef DEBUG_RDMon_AligTree
	printf(" *RDMon_AligTree: Evt %d Track %d: %d %d %d hits 0x%8x 0x%8x\n",
	       ptrEvt()->getEventNumberInRun(),
	       t.Id,nhAlig,nAlig[0],nAlig[1],iAlig[0],iAlig[1]);
#  endif
      if (nAlig[0]>TOpt::iCut[7] && nAlig[1]>TOpt::iCut[8] ||
	  nAlig[1]>TOpt::iCut[8] && nAlig[2]>TOpt::iCut[9] ||
	  nAlig[2]>TOpt::iCut[9] && nAlig[3]>TOpt::iCut[10]) {
	txAlig = t.Hfirst(3); tyAlig = t.Hfirst(4);
	AligTree->Fill();
#  ifdef DEBUG_RDMon_AligTree
	printf("OK\n");
#  endif
      }
    }   // End RDMon_AligTree
#endif
#ifdef RDMon_GEMWrBack
      if (nhGEMs_2>.75*nhAll_2) {
	GEMs_wrback = true;
	static int ievtp = -1;
	int ievt = CsEvent::Instance()->getEventNumberInRun();
	if (ievt!=ievtp) {
	  printf("Evt %8d RDMon_GEMWrBack  %d\n",ievt,t.Id);
	  ievtp = ievt;
	}
	else
	  printf("             RDMon_GEMWrBack %d\n",t.Id);
      }
      else if (nhGEMs_2>.3*nhAll_2) {
	printf("%d %d %d\n",t.Id,nhGEMs_2,nhAll_2);
      }
#endif
#ifdef RDMon_SM2Tree
      sm2_nhits = 0;
      for (int ipat = 0; ipat<SM2_NPATS; ipat++) sm2_pats[ipat] = 0;
      list<int>::const_iterator ih = t.lHPat().begin();
      while (ih!=t.lHPat().end()) {
	if (*ih>=0) {
	  THit& h = vecHit[*ih];
	  int ipl = h.IPlane;
	  int ipat = ipl/32;
	  if (ipat<SM2_NPATS) sm2_pats[ipat] |= (1<<(ipl%32));
	  sm2_nhits++;
	}
	ih++;
      }
      for (int i = 0; i<5; i++) sm2_hlx[i] = Hextr(i+1);
      sm2_chi2 = t.Chi2tot/(t.NHits-5);
      sm2_tree->Fill();
#endif
      //if (Hextr(0)>t.Hlast(0))
      //	t.Hlast = Hextr;                  // Reset "Hlast"
#ifdef DEBUG
      debug = 0;
#endif
      vsslat->Fill((double)vsslatz);

#ifdef DEBUG
      int extra2v = -1;
      if (extra2v!=-1) {
	if      (extra2v==-2) { // Extrapolate to ``vertex''
	  static double z = 0;
	  Hextr(0) = z; t.Hfirst.Extrapolate(Hextr); t.Hfirst = Hextr;
	}
	else if (extra2v==-3) {  // Delete track
	  listTrack.erase(it++); continue;
	}
	else if (extra2v==-4) {  // Split track @RICH
	  TTrack tx = *it; tx.Behead(112); // 112 = FI05Y
	  tx.FullKF(-1); tx.FullKF(1);
	  Hextr(0) = 750; tx.Hfirst.Extrapolate(Hextr); tx.Hfirst = Hextr;
	  listTrack.push_back(tx);
	  t.Shorten(113);
	  unsigned int beam_id = 0;
	  list<TTrack>::iterator jt = listTrack.begin(); 
	  if (beam_id) {
	    while (jt!=listTrack.end()) {  // Search for ID == "beam_id"
	      if ((*jt).Id==beam_id) break;
	      jt++;
	    }
	    if (jt==listTrack.end()) {
	      cout<<"DEBUG TEv::RDMonitor: No TTrack ID == "<<extra2v<<endl;
	      beam_id = 0;
	    }
	  }
	  if (beam_id) {
	    t.QuickKF(1,1); t.QuickKF(-1,1);
	    (*jt).Append(*it);
	    // Initial values
	    (*jt).Hfirst(5) = 1/160; (*jt).Hfirst(5,5) = 1.E-4;
	    (*jt).QuickKF(1,1); (*jt).QuickKF(-1,1);
	    (*jt).FullKF(1);  (*jt).FullKF(-1);
	    Hextr(0) = 750; (*jt).Hlast.Extrapolate(Hextr); (*jt).Hlast = Hextr;
	    listTrack.erase(it++); continue;
	  }
	  else { 
	    t.QuickKF(1,1); t.QuickKF(-1,1);
	    t.FullKF(1);    t.FullKF(-1);
	    Hextr(0) = 750; t.Hlast.Extrapolate(Hextr); t.Hlast = Hextr;
	  }
	}
	else if (extra2v>=0) {   // Append of ID == "extra2v"
	  list<TTrack>::iterator jt =listTrack.begin(); 
	  while (jt!=listTrack.end()) {  // Search for ID == "extra2v"
	    if ((int)(*jt).Id==extra2v) break;
	    jt++;
	  }
	  if (jt==listTrack.end()) {
	    cout<<"DEBUG TEv::RDMonitor: No TTrack ID == "<<extra2v<<endl;
	  }
	  else {
	    TTrack tx = *it; tx.Id = ++TTrack::TrackCounter;
	    tx.Append(*jt,"tmp");
	    int iter, ret;
	    for (iter = 0; iter<2; iter++) {
	      ret = iter==0 ? tx.QuickKF(-1,1) : tx.FullKF(-1); // backward
	      if(!ret) {
		CsErrLog::mes(elWarning,
			      "Backward Kalman fit failed => No Bridging");
		break;
	      }
	      ret = iter==0 ? tx.QuickKF( 1,1) : tx.FullKF( 1); // forward
	      if(!ret) {
		CsErrLog::mes(elWarning,
			      "Forward Kalman fit failed => No Bridging");
		break;
	      }
	    }
	    if (ret) {
	      t.Append(*jt);
	      for (iter = 0; iter<2; iter++) {
		ret = iter==0 ? t.QuickKF(-1,1) : t.FullKF(-1); // backward
		if(!ret) {
		  CsErrLog::mes(elWarning,
				"Backward Kalman fit failed => No Bridging");
		  break;
		}
		ret = iter==0 ? t.QuickKF( 1,1) : t.FullKF( 1); // forward
		if(!ret) {
		  CsErrLog::mes(elWarning,
				"Forward Kalman fit failed => No Bridging");
		  break;
		}
	      }
	    }
	  }
	}
      }
#endif
    }        // #################### End loop on tracks ####################

#ifdef RDMon_GEMWrBack
    if (GEMs_wrback && GEMs_stream>0) {
      CsEvent::Instance()->outputRawEventToStream(GEMs_stream);      
    }
#endif
#ifdef RDMon_DIGIT_GM
    // Noise on GEMs: single hit clusters not reconstructed
    for (igr = 0, idet = 0; igr<igr_mx; igr++)
      if ((1<<igr)&(TOpt::Hist[17]|    // Plane to monitor if inactive
		    TOpt::Hist[16])) { // Plane to monitor whatsoever
	for (ipl = setup.vIplFirst()[igr];
	     ipl<=setup.vIplLast()[igr]; ipl++, idet++) {
	  //                                              LOOP OVER PLANES
	  const TPlane &p = setup.vPlane(ipl);
	  const TDetect &d = setup.vDetect(p.IDetRef);
	  if (strncmp(d.Name.c_str(),"GM",2) &&
	      strncmp(d.Name.c_str(),"SI",2)) continue;// SELECT GEM/Si
	  int hist = gemhists[d.IDet];
	  vector<int>::const_iterator ihit;
	  for (ihit = p.vHitRef().begin(); ihit!=p.vHitRef().end(); ihit++) {
	    //                                            LOOP ON HITS
	    THit& h = vecHit[*ihit]; // ref to Thit vector element
	    if (h.sTrackID().size()) continue;         // NOT RECONSTRUCTED
	    list<CsDigit*> lDig = h.ptrCl->getDigitsList();
	    if (lDig.size()>1) continue;               // SINGLE HIT CLUSTER
	    CsDigit* Dig = lDig.front();
	    int wire = Dig->getAddress();
#  ifdef RDMon_GM_2D_CALIB
	    gemnoise[hist]->Fill((double)wire,Dig->getData()[2]);
#  else
	    gemnoise[hist]->Fill(Dig->getData()[2]);
#  endif
	  }
	}
      }
    
#endif
  }

#ifdef RDMon_ScifiSi_TIME
  if (TOpt::Hist[10] &&
      // ====================================================================
      // ********************** TIMING of SCIFI/SI HITS *********************
      // ====================================================================
      // Conditioned by trigger timing being precise enough:
      (evtTrig&stdTrigs) &&     // ``Standard'' triggers...
      !(evtTrig&highQ2Trig)) {  // ...excluding high-Q2 contribution

    static vector<TH2D*> hTPullsvsT, hTOffvsT, hTOffvsSiz;
    static vector<TH1D*> hTPulls, hTOffsets;
    static vector<int> pl2H_FI, pl2H_SI; static int ipl0;
    static bool first = true; static int nFIs; if (first) {
      //    ******************** HISTOGRAM BOOKING ********************
      first = false;
      if (setup.vIplFirst().size()<5)
	CsErrLog::mes(elFatal,
	  "\"TraF Hist[10]\" set while no scifi/Si zone is defined");
      CsHistograms::SetCurrentPath("/Traffic/RDmonitor");
      char hName[]   = "FI01X1_tOffvsSiz";
      //char hName[] = "FI01X1_tOffset";
      //char hName[] = "FI01X1_tPull";
      char hTitle[]   = "FI01X1 Time Offset w.r.t. Scifi Time  -  SIT#times!HighQ^{2}  ";
      //char hTitle[] = "FI01X1 Time Offset vs. |Track's Time|  ";
      //char hTitle[] = "FI01X1 Time Offset vs. Cluster Size  ";
      int ipl, hFI, hSI;
      for (ipl=ipl0 = setup.vIplFirst()[4], nFIs=hFI=hSI = 0;
	   ipl<=setup.vIplLast()[4]; ipl++) {
	const TPlane &p = setup.vPlane(ipl);
	const TDetect &d = setup.vDetect(p.IDetRef);
	if (d.IType!=21 && d.IType!=22) {
	  pl2H_FI.push_back(-1); pl2H_SI.push_back(-1); continue;
	}
	for (int i = 0; i<5; i++) { hName[i]=hTitle[i] = d.Name[i]; }
	if (d.IType==22) {
	  nFIs++;
	  pl2H_SI.push_back(-1); pl2H_FI.push_back(hFI++);
	  sprintf(hName+6,"_tPull"); sprintf(hTitle+6," T Pulls");
	  hTPulls.push_back(   new TH1D(hName,hTitle,100,-10,10));
	  sprintf(hName+6,"_tOffset"); sprintf(hTitle+6," Time");
	  hTOffsets.push_back( new TH1D(hName,hTitle,100,-10,10));
	  sprintf(hName+6,"_tOffvsSiz");
	  sprintf(hTitle+6," T vs. Cluster Size");
	  hTOffvsSiz.push_back(new TH2D(hName,hTitle,100,-10,10,8,.5,8.5));
	}
	else {
	  pl2H_FI.push_back(-1); pl2H_SI.push_back(hSI++);
	  sprintf(hName+6,"_tPull");
	  sprintf(hTitle+6," T-scifiT Pulls vs. |scifiT|");
	  hTPullsvsT.push_back(new TH2D(hName,hTitle,100,-10,10,5,0,10));
	  sprintf(hName+6,"_tOffset");
	  sprintf(hTitle+6," T-scifiT vs. |scifiT|");
	  hTOffvsT.push_back(  new TH2D(hName,hTitle,100,-10,10,5,0,10));
	}
      }
      CsHistograms::SetCurrentPath("/"); // ********** END of BOOKING **********
    }

    double siTSigScale = TOpt::dCut[79] ? TOpt::dCut[79] : 1;

    list<TTrack>::iterator it;
    for (it = listTrack.begin(); it!=listTrack.end(); it++) {
      TTrack &t = *it; if (!(t.Type&0x10)) continue;
      int nHits = t.NHits; if (nHits<TOpt::iPRpar[45]+2) continue;
      double chi2 = t.Chi2tot/(nHits-4); if (chi2>2) continue;

      // ********** LOOP ON GOOD (chi2<<, #hits>>) SCIFI/Si TRACKS **********

      int nScifis = 0; double s1 = 0, st = 0;
      list<int>::const_iterator ih = t.lHPat().begin();
      while (ih!=t.lHPat().end()) {
	if (*ih>=0) {
	  const THit &h = vecHit[*ih]; int ipl = h.IPlane;
	  const TPlane &p = setup.vPlane(ipl);
	  const TDetect &d = setup.vDetect(p.IDetRef);
	  if (d.IType==22) {
	    double hT = h.Time, w = 1/h.SigT/h.SigT; // Assuming time is always defined fo scifis..
	    nScifis++; s1 += w; st += hT*w;
	  }
	}
	ih++;
      }
      if (nFIs>=3 && nScifis<3 || nScifis<nFIs)      // ***** REQUIRE #SCIFIS >>
	continue;

      // ***** TAKE TRACK'S MEAN TIME from the SOLE SCIFIS
      // (The aim being mainly to check the Si's.)
      double tT = st/s1, tW2 = s1, tTa = fabs(tT);

      ih = t.lHPat().begin(); while (ih!=t.lHPat().end()) {
	if (*ih>=0) {
	  const THit &h = vecHit[*ih]; int jpl = h.IPlane-ipl0;
	  if (jpl>=(int)pl2H_FI.size()) {
	    printf("RDMon: Track #%d Plane #%d\n",t.Id,jpl+ipl0); break;
	  }
	  double hT = h.Time, sigt = h.SigT; // Assuming time is always defined for scifis and Sis...
	  sigt /= siTSigScale; double s = sqrt(sigt*sigt+1/tW2);
	  int hist = pl2H_FI[jpl]; if (hist>=0) {
	    //	    hTPulls[hist]->Fill(hT/h.SigT); hTOffsets[hist]->Fill(hT);
	    hTPulls[hist]->Fill((hT-tT)/s); hTOffsets[hist]->Fill(hT-tT);
	    hTOffvsSiz[hist]->Fill(hT,h.ptrCl->getDigitsList().size());
	  }
	  hist = pl2H_SI[jpl]; if (hist>=0) {
	    hTPullsvsT[hist]->Fill((hT-tT)/s,tTa);
	    hTOffvsT[hist]->Fill(hT-tT,tTa);
	  }
	}
	ih++;
      }
    }
  }  // End timing of scifi/si hits 
#endif
}
